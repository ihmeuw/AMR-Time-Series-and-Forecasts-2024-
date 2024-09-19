from __future__ import division

from builtins import range

import numpy as np
import pandas as pd
from sklearn import preprocessing

from db_queries import get_best_model_versions, get_covariate_estimates
from db_tools import ezfuncs

from cod_prep.claude import redistribution
from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import (
    add_age_metadata,
    add_cause_metadata,
    add_location_metadata,
    add_population,
    create_age_bins,
    get_current_location_hierarchy,
)
from cod_prep.utils import (
    drop_unmodeled_asc,
    drop_unmodeled_sex_causes,
    get_function_results,
    print_log_message,
    report_duplicates,
    report_if_merge_fail,
)

CONF = Configurator("standard")


def get_int_cause_model_covariates(int_cause):
    covariates_df = pd.read_csv("FILEPATH")
    assert len(covariates_df) == 1
    return covariates_df["covariates"].str.split(",").iloc[0]


def get_covariate_id(covariate):
    query = """SELECT * FROM shared.covariate WHERE covariate_name_short = '{}'""".format(covariate)
    df = ezfuncs.query(query, conn_def=CONF.get_database_setup("db"))
    assert len(df) == 1, f"{covariate} failed"
    return int(df["covariate_id"].iloc[0])


def get_covariate_id_cols(covariate_id):
    query = """SELECT * FROM shared.covariate WHERE covariate_id = {}""".format(covariate_id)
    df = ezfuncs.query(query, conn_def=CONF.get_database_setup("db"))
    assert len(df) == 1
    id_cols = ["location_id", "year_id"]
    if df["by_age"].iloc[0] == 1:
        id_cols += ["age_group_id"]
    if df["by_sex"].iloc[0] == 1:
        id_cols += ["sex_id"]
    return id_cols


def add_weight_to_covariate(cov_df, covariate_name_short, **cache_kwargs):
    if covariate_name_short is None:
        covariate_name_short = cov_df.loc[0, "covariate_name_short"]

    cov_to_weight_map = {
        "alcohol_lpc": "population",
        "SEV_metab_fpg": "population",
        "SEV_wash_sanitation": "population",
        "smok_prev": "population",
        "prev_obesity": "population",
        "idu_prop_agestandardized": "none",
        "idu_prop_byage": "population",
        "diabetes_fpg_byage": "population",
        "gbd_IP_env": "population",
        "smoking_prev": "population",
    }
    weight = cov_to_weight_map[covariate_name_short]

    if weight == "none":
        cov_df["weight"] = 1
    elif weight == "population":
        cov_df = add_population(
            cov_df,
            pop_run_id=CONF.get_id("pop_run"),
            location_set_id=CONF.get_id("location_set"),
            release_id=CONF.get_id("release"),
            **cache_kwargs,
        )
        report_if_merge_fail(
            cov_df,
            check_col="population",
            merge_cols=["age_group_id", "location_id", "year_id", "sex_id"],
        )
        cov_df.rename({"population": "weight"}, axis="columns", inplace=True)
    else:
        raise NotImplementedError

    return cov_df


def aggregate_covariate(
    cov_df, agg_col, covariate_name_short=None, agg_age_group_ids=None, id_cols=None, **cache_kwargs
):
    if covariate_name_short is None:
        covariate_name_short = cov_df.loc[0, "covariate_name_short"]
    if id_cols is None:
        covariate_id = get_covariate_id(covariate_name_short)
        id_cols = get_covariate_id_cols(covariate_id)

    assert agg_col in ["sex_id", "age_group_id"]

    cov_df = add_weight_to_covariate(cov_df, covariate_name_short, **cache_kwargs)

    cov_df["weighted_cov"] = cov_df["mean_value"] * cov_df["weight"]

    if agg_col == "age_group_id":
        cov_df = create_age_bins(cov_df, agg_age_group_ids, dropna=False)
    elif agg_col == "sex_id":
        assert set(cov_df.sex_id) == {1, 2}
        cov_df["sex_id"] = 3

    cov_df = cov_df.groupby(id_cols, as_index=False)[["weighted_cov", "weight"]].sum()

    cov_to_agg_type_map = {
        "alcohol_lpc": "average",
        "SEV_metab_fpg": "average",
        "SEV_wash_sanitation": "average",
        "smok_prev": "average",
        "prev_obesity": "average",
        "idu_prop_byage": "average",
        "idu_prop_agestandardized": "average",
        "diabetes_fpg_byage": "average",
        "gbd_IP_env": "average",
        "smoking_prev": "average",
    }
    agg_type = cov_to_agg_type_map[covariate_name_short]

    if agg_type == "average":
        cov_df.eval("agg_cov = weighted_cov / weight", inplace=True)
    elif agg_type == "sum":
        cov_df["agg_cov"] = cov_df["weighted_cov"]
    else:
        raise NotImplementedError

    cov_df.drop(["weighted_cov", "weight"], axis="columns", inplace=True)
    cov_df.rename({"agg_cov": "mean_value"}, axis="columns", inplace=True)
    assert set(cov_df.columns) == set(id_cols + ["mean_value"])
    assert cov_df.notnull().values.all()

    return cov_df


def merge_covariate(df, covariate_name_short, scale=False, **get_cov_kwargs):
    covariate_id = get_covariate_id(covariate_name_short)
    id_cols = get_covariate_id_cols(covariate_id)
    get_cov_kwargs.update(
        {"location_id": list(df.location_id.unique()), "year_id": list(df.year_id.unique())}
    )
    cov_df = get_cov(covariate_id=covariate_id, **get_cov_kwargs)[id_cols + ["mean_value"]]
    cache_kwargs = {"force_rerun": False, "block_rerun": True}

    if "sex_id" in id_cols:
        if not set(df.sex_id).issubset(set(cov_df.sex_id)):
            assert set(df.sex_id) == {3}
            assert set(cov_df.sex_id) == {1, 2}
            cov_df = aggregate_covariate(
                cov_df,
                "sex_id",
                covariate_name_short=covariate_name_short,
                id_cols=id_cols,
                **cache_kwargs,
            )

    if "age_group_id" in id_cols:
        df_age_group_ids = df.age_group_id.unique().tolist()
        cov_age_group_ids = cov_df.age_group_id.unique().tolist()

        if not set(df_age_group_ids).issubset(set(cov_age_group_ids)):
            print_log_message("Aggregating covariates to match incoming dataframe.")
            cov_df = add_age_metadata(
                cov_df, ["age_group_days_start", "age_group_days_end"], **cache_kwargs
            )
            df_ag = add_age_metadata(
                df.copy(), ["age_group_days_start", "age_group_days_end"], **cache_kwargs
            )
            too_young = cov_df["age_group_days_end"] <= df_ag.age_group_days_start.min()
            too_old = cov_df["age_group_days_start"] >= df_ag.age_group_days_end.max()
            cov_df = cov_df[~(too_young | too_old)]
            cov_df = cov_df.drop(["age_group_days_start", "age_group_days_end"], axis=1)
            cov_df = aggregate_covariate(
                cov_df,
                "age_group_id",
                covariate_name_short=covariate_name_short,
                agg_age_group_ids=df_age_group_ids,
                id_cols=id_cols,
                **cache_kwargs,
            )

    cov_df = cov_df.rename(columns={"mean_value": covariate_name_short})
    report_duplicates(cov_df, id_cols)
    print_log_message("RUNNING MERGE COVARIATES")
    if covariate_name_short == "haqi":
        cov_df.eval(f"{covariate_name_short} = {covariate_name_short} / 100", inplace=True)
    if (covariate_name_short == "LDI_pc") or ("ldi_pc" in covariate_name_short):
        print_log_message("using natural log of LDI per capita")
        cov_df[covariate_name_short] = np.log(cov_df[covariate_name_short])

    if scale:
        scaler = preprocessing.MinMaxScaler()
        cov_df[[covariate_name_short]] = scaler.fit_transform(cov_df[[covariate_name_short]])

    df = df.merge(cov_df, on=id_cols, how="left")
    report_if_merge_fail(df, covariate_name_short, id_cols)
    return df


def drop_unmodeled_asc_aggregate_ages(df, cause_meta_df, age_meta_df, detail_ages):
    simple_age_map = age_meta_df.loc[lambda d: d["age_group_id"].isin(detail_ages), :].set_index(
        "simple_age", verify_integrity=True
    )
    df = add_cause_metadata(
        df,
        add_cols=["cause_start", "cause_end", "male", "female", "yld_only"],
        cause_meta_df=cause_meta_df.assign(
            cause_start=lambda d: d["yll_age_start"].map(simple_age_map["age_group_years_start"]),
            cause_end=lambda d: d["yll_age_end"].map(simple_age_map["age_group_years_end"]),
        ),
    )
    df = add_age_metadata(
        df, add_cols=["age_group_years_start", "age_group_years_end"], age_meta_df=age_meta_df
    )
    df = drop_unmodeled_sex_causes(df)

    unmodeled_age_old = df["age_group_years_start"] >= df["cause_end"]
    unmodeled_age_young = df["age_group_years_end"] <= df["cause_start"]
    df = df[~(unmodeled_age_old | unmodeled_age_young | (df["yld_only"] == 1))]
    df.drop(
        [
            "cause_start",
            "cause_end",
            "male",
            "female",
            "yld_only",
            "age_group_years_start",
            "age_group_years_end",
        ],
        axis=1,
        inplace=True,
    )
    return df


def enforce_redistribution_restrictions(int_cause, cause_meta_df, age_meta_df, template):
    ignore_acauses = []

    target_restriction_map = (
        pd.read_csv("FILEPATH" if int_cause in {"x59", "y34"} else "FILEPATH")
        .loc[lambda d: ~d["acause"].isin(ignore_acauses), :]
        .pipe(
            redistribution.process_redistribution_target_restrictions,
            cause_meta_df=cause_meta_df,
            allow_distribution_onto_yld_only=False,
            allow_logic_overrides_to_override_cause_hierarchy=False,
        )
        .set_index("cause_id")["restrictions"]
    )
    loc_meta_df = get_current_location_hierarchy(
        CONF.get_id("location_set_version"),
        force_rerun=False,
        block_rerun=True,
    ).assign(
        country_id=lambda d: d["path_to_top_parent"]
        .str.split(",")
        .map(lambda x: x[3] if len(x) > 3 else None)
        .astype(float),
    )
    age_meta_df = age_meta_df.rename(columns={"simple_age": "age"})
    return pd.concat(
        group.assign(evaluated=group.eval(target_restriction_map.get(cause_id, "True")))
        .query("evaluated")
        .loc[:, template.columns]
        for cause_id, group in (
            template.pipe(add_location_metadata, "country_id", location_meta_df=loc_meta_df)
            .pipe(add_age_metadata, "age", age_meta_df=age_meta_df)
            .groupby("cause_id")
        )
    )


def get_cov(
    covariate_id=None,
    model_version_id=None,
    location_set_id=None,
    release_id=None,
    location_id="all",
    year_id="all",
    **cache_kwargs,
):
    if location_set_id is None:
        location_set_id = CONF.get_id("location_set")
    if release_id is None:
        release_id = CONF.get_id("release")
    if model_version_id is None:
        model_version_id = get_best_model_versions(
            entity="covariate",
            ids=covariate_id,
            release_id=release_id,
            status="best",
        ).loc[0, "model_version_id"]
    cache_name = f"cov_{covariate_id}_mvid_{model_version_id}_lsid_{location_set_id}"
    function = get_covariate_estimates
    args = [covariate_id]
    kwargs = {
        "location_set_id": location_set_id,
        "release_id": release_id,
        "model_version_id": model_version_id,
        "location_id": location_id,
        "year_id": year_id,
    }
    df = get_function_results(function, args, kwargs, cache_name, **cache_kwargs)
    if isinstance(df, pd.DataFrame):
        assert (
            covariate_id == df.loc[0, "covariate_id"]
        ), f"Covariate {covariate_id} and model version {model_version_id} do not match."
    return df
