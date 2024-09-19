import logging
import os
import sys
from builtins import range, str

import numpy as np
import pandas as pd

from db_queries import get_outputs
from get_draws.api import get_draws

from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import (
    create_age_bins,
    get_ages,
    get_cause_map,
    get_current_cause_hierarchy,
    get_current_location_hierarchy,
)
from cod_prep.downloaders.causes import get_all_related_causes
from cod_prep.utils import print_log_message
from mcod_prep.mcod_mapping import MCoDMapper
from mcod_prep.utils.causes import get_all_related_syndromes, get_int_cause_hierarchy
from mcod_prep.utils.mcause_io import McauseResult, setup_logger

CONF = Configurator("standard")
DEM_COLS = ["cause_id", "location_id", "sex_id", "year_id", "age_group_id"]
DRAW_COLS = ["draw_" + str(x) for x in range(0, CONF.get_id("draw_num"))]
EST_COLS = DRAW_COLS + ["point_estimate"]
BLOCK_RERUN = {"force_rerun": False, "block_rerun": True}


def convert_int_cols(df):
    df[DEM_COLS] = df[DEM_COLS].apply(pd.to_numeric, downcast="integer")
    return df


def get_location_list():
    lhh = get_current_location_hierarchy(
        location_set_version_id=CONF.get_id("location_set_version"), **BLOCK_RERUN
    ).query("most_detailed == 1")
    return list(lhh["location_id"].unique())


def get_codcorrect_deaths(
    cause_id,
    year_id,
    age_group_id,
    location_id,
    sex_id,
    logger,
    num_workers=None,
    gbd_round_id=None,
    release_id=None,
    version_id=None,
    downsample=True,
    n_draws=CONF.get_id("draw_num"),
):
    if gbd_round_id is None:
        gbd_round_id = CONF.get_id("gbd_round")
    if release_id is None:
        release_id = CONF.get_id("release")
    if version_id is None:
        version_id = CONF.get_id("codcorrect_version")
    get_draws_kwargs = {
        "gbd_id_type": "cause_id",
        "gbd_id": int(cause_id),
        "source": "codcorrect",
        "metric_id": 1,
        "measure_id": 1,
        "location_id": location_id,
        "age_group_id": age_group_id,
        "year_id": int(year_id),
        "version_id": int(version_id),
        "num_workers": num_workers,
        "release_id": int(release_id),
        "sex_id": sex_id,
        "downsample": downsample,
        "n_draws": n_draws,
    }
    get_outputs_kwargs = {
        "topic": "cause",
        "cause_id": int(cause_id),
        "metric_id": 1,
        "measure_id": 1,
        "location_id": location_id,
        "age_group_id": age_group_id,
        "year_id": int(year_id),
        "release_id": int(release_id),
        "compare_version_id": CONF.get_id("compare_version"),
        "sex_id": sex_id,
    }
    logger.debug("Reading draws of deaths from CoDCorrect")
    df_draws = get_draws(**get_draws_kwargs)

    logger.debug("Reading output of deaths from CoDCorrect")
    df_point = get_outputs(**get_outputs_kwargs)
    df_point = df_point.rename(columns={"val": "point_estimate"})

    logger.debug("Merging point estimate onto draws")
    df = df_draws.merge(
        df_point[DEM_COLS + ["point_estimate"]], how="left", on=DEM_COLS, validate="one_to_one"
    )

    return df


def get_deaths_draws(
    cause_id,
    year_id,
    ages,
    locs,
    sexes,
    int_cause,
    end_product,
    logger,
    gbd_round_id=None,
    release_id=None,
    version_id=None,
    get_all_ages=True,
    downsample=True,
    n_draws=CONF.get_id("draw_num"),
):
    if gbd_round_id is None:
        gbd_round_id = CONF.get_id("gbd_round")
    if release_id is None:
        release_id = CONF.get_id("release")
    if version_id is None:
        version_id = CONF.get_id("codcorrect_version")
    is_syndrome = int_cause in MCoDMapper.infectious_syndromes
    get_codcorrect_kwargs = {
        "cause_id": cause_id,
        "location_id": locs,
        "age_group_id": ages,
        "year_id": year_id,
        "gbd_round_id": gbd_round_id,
        "version_id": version_id,
        "release_id": release_id,
        "sex_id": sexes,
        "downsample": downsample,
        "n_draws": n_draws,
        "logger": logger,
    }
    if end_product in ["incidence", "attributable_burden"]:
        mort_result = McauseResult(
            int_cause=int_cause,
            end_product="mortality",
            process="calculate_counts",
            year_id=year_id,
            cause_id=cause_id,
            status="latest",
            conf=CONF,
        )
        logger.debug(f"Reading {int_cause} deaths from\n {mort_result.parent_model_dir}")
        df = mort_result.read_results()
        if end_product == "attributable_burden":
            logger.debug("Reading draws of deaths from CoDCorrect")
            cc_df = get_codcorrect_deaths(**get_codcorrect_kwargs)
            df = df.merge(cc_df, how="left", on=DEM_COLS)
            logger.debug(f"Calculating non {int_cause}-related deaths")
            for draw in DRAW_COLS:
                df[draw] = df[f"{draw}_y"] - df[f"{draw}_x"]
            assert (df[DRAW_COLS] >= 0).values.all()
    elif end_product == "mortality" and is_syndrome:
        sepsis_result = McauseResult(
            int_cause="sepsis",
            end_product="mortality",
            process="cause_aggregation",
            year_id=year_id,
            cause_id=cause_id,
            description=CONF.get_id("sepsis_model"),
            conf=CONF,
        )
        logger.debug(f"Reading sepsis deaths from\n {sepsis_result.parent_model_dir}")
        df = sepsis_result.read_results()
    elif end_product in ["mortality", "rdp"] and not is_syndrome:
        logger.debug("Reading draws of deaths from CoDCorrect")
        df = get_codcorrect_deaths(**get_codcorrect_kwargs)
    missing_ages = set(ages) - set(df.age_group_id)
    if get_all_ages:
        assert missing_ages == set(), f"Missing ages: {missing_ages}"
    else:
        if missing_ages != set():
            logger.debug(f"Missing the following age groups from the denominator: {missing_ages}")
    assert df[DEM_COLS + EST_COLS].notnull().values.all(), "Error pulling draws of deaths"
    return df.groupby(DEM_COLS, as_index=False)[EST_COLS].sum()


def convert_prop_to_count(prop_df, deaths_df, end_product, logger):
    logger.debug("Calculating counts")
    if end_product in ["mortality", "rdp"]:
        df = prop_df.merge(deaths_df, how="left", on=DEM_COLS)
    elif end_product in ["incidence", "attributable_burden"]:
        prop_df[DRAW_COLS] = np.reciprocal(prop_df[DRAW_COLS])
        df = prop_df.merge(deaths_df, how="inner", on=DEM_COLS)
    for draw in EST_COLS:
        df[draw] = df[f"{draw}_x"] * df[f"{draw}_y"]
    df = df[DEM_COLS + EST_COLS]
    assert df.notnull().values.all(), "Error calculating deaths"
    return df


def squeeze(description, int_cause, end_product, year_id, cause_id, prop_df, logger):
    chh = get_current_cause_hierarchy(
        cause_set_version_id=CONF.get_id("computation_cause_set_version"),
        cause_set_id=CONF.get_id("computation_cause_set"),
        force_rerun=False,
        block_rerun=True,
    )

    scl = pd.read_csv("FILEPATH")

    if int_cause == "sepsis":
        all_sepsis_cause_ids = scl.loc[scl["100_sepsis"] == 1, "cause_id"].unique()
        all_sepsis_causes = []
        for cause in all_sepsis_cause_ids:
            sepsis_causes = get_all_related_causes(cause, chh)
            all_sepsis_causes = all_sepsis_causes + [
                x for x in sepsis_causes if x not in all_sepsis_causes
            ]

        if cause_id in all_sepsis_causes:
            prop_df[EST_COLS] = 1
    else:
        version_df = pd.read_csv("FILEPATH").query("int_cause != @int_cause")

        prop_df["int_cause"] = int_cause
        prediction_dfs = [prop_df]

        syn_cause_map = pd.melt(
            scl,
            id_vars=["cause_id", "level", "acause"],
            value_vars=[x for x in scl.columns if "L2_" in x or "L1_" in x or "L1.5_" in x],
            var_name="infectious_syndrome",
            value_name="modeled",
        )
        modeled_syndromes = (
            syn_cause_map.loc[
                (syn_cause_map["cause_id"] == cause_id)
                & (syn_cause_map["modeled"] == 1)
                & (syn_cause_map["infectious_syndrome"] != int_cause),
                "infectious_syndrome",
            ]
            .unique()
            .tolist()
        )

        for syndrome in modeled_syndromes:
            logger.debug(f"Appending infectious syndrome: {syndrome}")
            predict_description = "".join(
                version_df.loc[version_df["int_cause"] == syndrome, "model_description"].tolist()
            )

            predict_df = McauseResult(
                int_cause=syndrome,
                end_product=end_product,
                process="predict",
                year_id=year_id,
                cause_id=cause_id,
                description=predict_description,
                conf=CONF,
            ).read_results()

            if "level_2" in predict_df.columns:
                predict_df.rename(columns={"level_2": "cause_id"}, inplace=True)
                convert_int_cols(predict_df)

            predict_df["int_cause"] = syndrome

            prediction_dfs.append(predict_df)

        prediction_dfs = pd.concat(prediction_dfs, ignore_index=True, sort=True)

        prediction_dfs = prediction_dfs[
            [x for x in list(prediction_dfs) if ("_id" in x)]
            + ["level_1", "int_cause", "point_estimate"]
            + [x for x in list(prediction_dfs) if ("draw" in x)]
        ]

        assert set(prediction_dfs["int_cause"].unique().tolist()) == set(
            modeled_syndromes + [int_cause]
        )

        predict_cols = [x for x in list(prediction_dfs) if ("id" in x)] + ["level_1"]

        prediction_dfs[EST_COLS] = pd.DataFrame(
            prediction_dfs[EST_COLS].to_numpy()
            / prediction_dfs.groupby(predict_cols)[EST_COLS].transform(sum).to_numpy(),
            index=prediction_dfs.index,
        )

        assert (
            np.unique(
                np.isclose(
                    prediction_dfs.groupby(predict_cols)[EST_COLS].transform(sum).to_numpy(), 1
                )
            )
            == True
        )

        prediction_dfs = prediction_dfs.loc[prediction_dfs["int_cause"] == int_cause]
        prediction_dfs = prediction_dfs.drop(["int_cause"], axis=1)

        prop_df = prediction_dfs

    logger.debug("Saving squeeze")
    McauseResult(
        int_cause=int_cause,
        end_product=end_product,
        process="squeeze",
        year_id=year_id,
        cause_id=cause_id,
        description=description,
        conf=CONF,
    ).write_results(prop_df)

    return prop_df


def main_calculate_counts(
    description, cause_id, year_id, end_product, int_cause, no_squeeze, logger=None
):
    if logger == None:
        logger = setup_logger(f"calculate_counts_{year_id}_{cause_id}")

    prop_df = McauseResult(
        int_cause=int_cause,
        end_product=end_product,
        process="predict",
        year_id=year_id,
        cause_id=cause_id,
        description=description,
        conf=CONF,
    ).read_results()
    if "level_2" in prop_df.columns:
        prop_df.rename(columns={"level_2": "cause_id"}, inplace=True)
    convert_int_cols(prop_df)
    if (
        not no_squeeze
        and int_cause in MCoDMapper.infectious_syndromes + ["sepsis"]
        and end_product == "mortality"
    ):
        prop_df = squeeze(description, int_cause, end_product, year_id, cause_id, prop_df, logger)

    ages = list(prop_df.age_group_id.unique())
    locs = list(prop_df.location_id.unique())
    sexes = list(prop_df.sex_id.unique())
    draws_df = get_deaths_draws(
        cause_id, year_id, ages, locs, sexes, int_cause, end_product, logger
    )
    df = convert_prop_to_count(prop_df, draws_df, end_product, logger)
    logger.debug("Saving output")
    McauseResult(
        int_cause=int_cause,
        end_product=end_product,
        process="calculate_counts",
        year_id=year_id,
        cause_id=cause_id,
        description=description,
        conf=CONF,
    ).write_results(df)


if __name__ == "__main__":
    description = str(sys.argv[1])
    int_cause = str(sys.argv[2])
    end_product = str(sys.argv[3])
    no_squeeze = str(sys.argv[4])

    assert no_squeeze in ["True", "False"]
    no_squeeze = no_squeeze == "True"

    task_id = os.getenv("SLURM_ARRAY_TASK_ID")
    if task_id:
        print_log_message(f"Running as array job, task_id: {task_id}")
        base_out_dir = McauseResult(
            int_cause=int_cause,
            end_product=end_product,
            process="run_model",
            description=description,
            conf=CONF,
        ).parent_model_dir
        task_row = pd.read_parquet("FILEPATH").iloc[int(task_id) - 1]
        cause_id = int(task_row["cause_id"])
        year_id = int(task_row["year_id"])
    else:
        cause_id = int(sys.argv[5])
        year_id = int(sys.argv[6])
    print_log_message(f"Running year {year_id}, cause {cause_id}")
    print_log_message(f"No Squeeze: {no_squeeze}")
    main_calculate_counts(description, cause_id, year_id, end_product, int_cause, no_squeeze)
