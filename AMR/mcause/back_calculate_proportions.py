import getpass
import logging
import os
import sys
from builtins import range, str

import numpy as np
import pandas as pd
from scipy import stats

from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import get_cod_ages
from cod_prep.utils import print_log_message, report_if_merge_fail
from mcod_prep.age_loc_aggregation import create_age_standardized_rates
from mcod_prep.calculate_counts import get_deaths_draws
from mcod_prep.utils.mcause_io import McauseResult, makedirs_safely, setup_logger

CONF = Configurator("standard")
DEM_COLS = ["cause_id", "location_id", "sex_id", "year_id", "age_group_id"]
DRAW_COLS = ["draw_" + str(x) for x in range(0, CONF.get_id("draw_num"))]
EST_COLS = DRAW_COLS + ["point_estimate"]
BLOCK_RERUN = {"force_rerun": False, "block_rerun": True}


def get_summary_stats(df, cols, suffix):
    df = df.copy()
    df[f"mean_{suffix}"] = np.mean(df[cols].values.flatten())
    df[f"median_{suffix}"] = np.median(df[cols].values.flatten())
    df[f"min_{suffix}"] = np.min(df[cols].values.flatten())
    df[f"max_{suffix}"] = np.max(df[cols].values.flatten())
    return df


def calculate_prop(df, draw, parent_model_dir, year_id, cause_id, logger):
    over_one = df[f"{draw}_x"] > df[f"{draw}_y"]
    outside_tol = ~np.isclose(df[f"{draw}_x"], df[f"{draw}_y"], rtol=0.01, atol=10)
    bad_rows = over_one & outside_tol
    if bad_rows.any():
        bad = df[bad_rows]
        x_cols = [x for x in bad.columns if "_x" in x]
        y_cols = [x for x in bad.columns if "_y" in x]
        bad = get_summary_stats(bad, x_cols, "x")
        bad = get_summary_stats(bad, y_cols, "y")
        bad = bad[
            DEM_COLS
            + [
                f"{draw}_x",
                f"{draw}_y",
                "point_estimate_x",
                "point_estimate_y",
                "mean_x",
                "mean_y",
                "median_x",
                "median_y",
                "min_x",
                "min_y",
                "max_x",
                "max_y",
            ]
        ]

        bad.to_csv("FILEPATH", index=False)

    if draw == "point_estimate":
        assert not bad_rows.any(), (
            f"Found the following bad rows in {draw} with proportion "
            f"greater than 1: \n {df[bad_rows]}"
        )

    return np.clip((df[f"{draw}_x"] / df[f"{draw}_y"]).fillna(0), None, 1)


def convert_count_to_prop(counts_df, deaths_df, parent_model_dir, year_id, cause_id, logger):
    logger.debug("Calculating proportions")
    df = counts_df.merge(deaths_df, how="left", on=DEM_COLS)
    report_if_merge_fail(df, EST_COLS[0] + "_y", DEM_COLS)
    for draw in EST_COLS:
        df[draw] = calculate_prop(df, draw, parent_model_dir, year_id, cause_id, logger)
    df = df[DEM_COLS + EST_COLS]
    logger.debug("Running checks")
    assert (
        df.notnull().values.all()
    ), f"Error calculating proportions \n{df.columns[df.isnull().any()]}"
    assert (df[EST_COLS] <= 1).all().all(), f"Proportion should not exceed 1"
    return df


def main_props(
    description, cause_id, year_id, end_product, int_cause, parent_model_dir, logger=None
):
    if logger == None:
        logger = setup_logger(f"back_calculate_props_{int_cause}_{cause_id}_{year_id}")

    counts_df = McauseResult(
        int_cause=int_cause,
        end_product=end_product,
        process="cause_aggregation",
        year_id=year_id,
        cause_id=cause_id,
        description=description,
        conf=CONF,
    ).read_results()
    counts_df[DEM_COLS] = counts_df[DEM_COLS].astype(int)
    counts_df = counts_df.query("age_group_id != 27")
    ages = list(counts_df.age_group_id.unique())
    locs = list(counts_df.location_id.unique())
    sexes = list(counts_df.sex_id.unique())
    deaths_df = get_deaths_draws(
        cause_id,
        year_id,
        ages,
        locs,
        sexes,
        int_cause,
        end_product,
        logger=logger,
        get_all_ages=False,
    )
    df = convert_count_to_prop(counts_df, deaths_df, parent_model_dir, year_id, cause_id, logger)
    logger.debug("Creating age standardized rates")
    age_std_df = create_age_standardized_rates(df, EST_COLS, DEM_COLS)
    df = pd.concat([age_std_df, df], ignore_index=True, sort=True)
    logger.debug("Saving output")
    McauseResult(
        int_cause=int_cause,
        end_product=end_product,
        process="calculate_props",
        year_id=year_id,
        cause_id=cause_id,
        description=description,
        conf=CONF,
    ).write_results(df)


if __name__ == "__main__":
    description = str(sys.argv[1])
    int_cause = str(sys.argv[2])
    end_product = str(sys.argv[3])
    parent_model_dir = str(sys.argv[4])
    assert end_product != "rdp", "Not used for redistribution props"
    task_id = os.getenv("SLURM_ARRAY_TASK_ID")
    if task_id:
        print(f"Running as array job, task_id: {task_id}")
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
    main_props(description, cause_id, year_id, end_product, int_cause, parent_model_dir)
