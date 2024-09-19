import argparse
import os
import sys
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd

from get_draws.api import get_draws

from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import (
    add_population,
    get_cod_ages,
    get_current_cause_hierarchy,
    pretty_print,
)
from cod_prep.utils import print_log_message, report_if_merge_fail
from mcod_prep.utils.mcause_io import McauseResult, setup_logger

CONF = Configurator("standard")
DEM_COLS = ["cause_id", "location_id", "sex_id", "year_id", "age_group_id"]
DRAW_COLS = ["draw_" + str(x) for x in range(0, CONF.get_id("draw_num"))]
EST_COLS = DRAW_COLS + ["point_estimate"]


def convert_int_cols(df):
    df[DEM_COLS] = df[DEM_COLS].apply(pd.to_numeric, downcast="integer")
    return df


def get_isclose(x):
    isclose = np.isclose(x[0], x[1], rtol=0.00001)
    return isclose


def summarize_draws(df, draw_cols=DRAW_COLS, prefix=""):
    constant = 0.00000000000001

    df[draw_cols] = np.log(df[draw_cols] + constant)

    df[prefix + "log_mean"] = np.log(df["point_estimate"] + constant)

    df["log_std_dev"] = np.std(df[draw_cols], axis=1, ddof=1)

    df[prefix + "upper"] = np.exp(df[prefix + "log_mean"] + (1.96 * df["log_std_dev"])) - constant
    df[prefix + "lower"] = np.exp(df[prefix + "log_mean"] - (1.96 * df["log_std_dev"])) - constant

    df[prefix + "mean"] = np.exp(df[prefix + "log_mean"]) - constant

    if prefix != "rate_":
        assert (
            df[prefix + "lower"] > -0.00000001
        ).all(), "Not all lower bound values are greater than 0."
    else:
        assert (
            df.loc[df["age_group_id"] != 27][prefix + "lower"] > -0.00000001
        ).all(), "Not all lower bound values are greater than 0."

    df["is_close_mean"] = df[["point_estimate", prefix + "mean"]].apply(
        lambda x: get_isclose(x), axis=1
    )

    if prefix == "rate_":
        df.loc[df["age_group_id"] == 27, "is_close_mean"] = True
        df.loc[df["age_group_id"] == 27, "is_close_std"] = True

    condition_failed_mean = df[df["is_close_mean"] != True]

    if not condition_failed_mean.empty:
        print("Assertion failed for the following rows (mean):")
        print(condition_failed_mean)

    if condition_failed_mean.empty:
        print("Assertion passed. All values in 'mean_is_close' are True.")

    drop_cols = ["point_estimate"] + draw_cols

    drop_cols += [prefix + "log_mean", "log_std_dev"]
    df = df.drop(columns=drop_cols, errors="ignore")

    return df


def compile_metrics(cause, year, description, end_product, int_cause, logger):
    logger.debug(f"Reading in cause {cause} count files")
    df = convert_int_cols(
        McauseResult(
            int_cause=int_cause,
            end_product=end_product,
            process="cause_aggregation",
            description=description,
            year_id=year,
            cause_id=cause,
            conf=CONF,
        ).read_results()
    )

    if end_product == "mortality":
        logger.debug(f"Calculating cause {cause} specific chain fraction")
        denom_df = convert_int_cols(
            McauseResult(
                int_cause=int_cause,
                end_product=end_product,
                process="cause_aggregation",
                description=description,
                year_id=year,
                cause_id=294,
                conf=CONF,
            ).read_results()
        ).drop("cause_id", axis=1)
        cscf_df = df.merge(denom_df, on=[x for x in DEM_COLS if x != "cause_id"], how="left")
        cscf_df[DRAW_COLS] = cscf_df.filter(regex="draw_.*_x").rename(
            columns=lambda c: c[:-2]
        ) / cscf_df.filter(regex="draw_.*_y").rename(columns=lambda c: c[:-2])
        cscf_df["point_estimate"] = cscf_df["point_estimate_x"] / cscf_df["point_estimate_y"]
        cscf_df = summarize_draws(cscf_df, prefix="cscf_")[
            DEM_COLS + ["cscf_mean", "cscf_upper", "cscf_lower"]
        ]

    logger.debug(f"Reading in cause {cause} proportion files")
    cf_df = McauseResult(
        int_cause=int_cause,
        end_product=end_product,
        process="calculate_props",
        description=description,
        year_id=year,
        cause_id=cause,
        conf=CONF,
    ).read_results()
    cf_df = convert_int_cols(cf_df)
    prefix = "cf_"
    if end_product == "incidence":
        prefix = "cfr_"
    cf_df = summarize_draws(cf_df, prefix=prefix)[DEM_COLS + ["cf_mean", "cf_upper", "cf_lower"]]

    logger.debug(f"Preparing rates for cause {cause}")
    rate_df = df.copy()
    rate_df = convert_int_cols(rate_df)
    prefix = "rate_"

    ages = rate_df["age_group_id"].unique().tolist()
    locations = rate_df["location_id"].unique().tolist()
    years = rate_df["year_id"].unique().tolist()
    sexes = rate_df["sex_id"].unique().tolist()

    cache_kwargs = {"force_rerun": False, "block_rerun": True}
    pop_run_id = CONF.get_id("pop_run")
    rate_df = add_population(rate_df, pop_run_id=pop_run_id, **cache_kwargs)
    report_if_merge_fail(rate_df.loc[rate_df["age_group_id"] != 27], "population", DEM_COLS)

    for col in EST_COLS:
        rate_df[col] = rate_df[col] / rate_df["population"]
    rate_df = summarize_draws(rate_df, prefix=prefix)[
        DEM_COLS + ["rate_mean", "rate_upper", "rate_lower"]
    ]

    logger.debug(f"Merging {cause} metrics together")
    df = summarize_draws(df)[DEM_COLS + ["mean", "upper", "lower"]]
    merge_df = (
        df.merge(cf_df, on=DEM_COLS, how="outer", validate="one_to_one")
        .merge(cscf_df, on=DEM_COLS, how="outer", validate="one_to_one")
        .merge(rate_df, on=DEM_COLS, how="outer", validate="one_to_one")
    )

    report_if_merge_fail(merge_df, "mean", DEM_COLS)
    report_if_merge_fail(merge_df, "cf_mean", DEM_COLS)
    report_if_merge_fail(merge_df, "cscf_mean", DEM_COLS)

    assert (
        merge_df[
            [x for x in merge_df.columns if x not in ["rate_mean", "rate_upper", "rate_lower"]]
        ]
        .notnull()
        .values.all()
    )

    condition = merge_df["age_group_id"] == 27
    merge_df.loc[condition, "rate_mean"] = merge_df["mean"]
    merge_df.loc[condition, "rate_upper"] = merge_df["upper"]
    merge_df.loc[condition, "rate_lower"] = merge_df["lower"]

    merge_df.loc[
        condition,
        [
            "mean",
            "upper",
            "lower",
            "cf_mean",
            "cf_upper",
            "cf_lower",
            "cscf_mean",
            "cscf_upper",
            "cscf_lower",
        ],
    ] = np.nan

    logger.debug(f"Job done: {cause}")
    return merge_df


def main_compile(description, end_product, int_cause, year_id, cause_id, logger=None):
    if logger == None:
        logger = setup_logger(f"compile_{int_cause}_{cause_id}_{year_id}")

    df = compile_metrics(cause_id, year_id, description, end_product, int_cause, logger)
    McauseResult(
        int_cause=int_cause,
        end_product=end_product,
        process="compile",
        description=description,
        year_id=year_id,
        cause_id=cause_id,
        conf=CONF,
    ).write_results(df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile summary files")
    parser.add_argument("description", type=str)
    parser.add_argument("end_product", type=str)
    parser.add_argument("int_cause", type=str)
    parser.add_argument("cause_id", type=int)
    parser.add_argument("year_id", type=int)
    args = parser.parse_args()
    main_compile(**vars(args))
