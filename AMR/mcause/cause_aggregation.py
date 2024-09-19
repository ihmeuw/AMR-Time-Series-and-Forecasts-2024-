import argparse
import logging
import os
import sys
from builtins import range, str

import pandas as pd

from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import get_all_child_cause_ids, get_current_cause_hierarchy
from cod_prep.utils import print_log_message, report_duplicates
from mcod_prep.age_loc_aggregation import create_age_standardized_rates
from mcod_prep.utils.mcause_io import McauseResult, setup_logger

CONF = Configurator("standard")
DEM_COLS = ["cause_id", "location_id", "sex_id", "year_id", "age_group_id"]
DRAW_COLS = ["draw_" + str(x) for x in range(0, CONF.get_id("draw_num"))]
EST_COLS = DRAW_COLS + ["point_estimate"]
CACHE_KWARGS = {"force_rerun": False, "block_rerun": True}


def main_cause_agg(description, year_id, int_cause, end_product, parent_cause_id, logger=None):
    if logger == None:
        logger = setup_logger(f"cause_aggregation_{int_cause}_{parent_cause_id}_{year_id}")

    logger.debug("Reading in age/sex/location aggregated files")

    chh = get_current_cause_hierarchy(
        cause_set_version_id=CONF.get_id("cause_set_version"),
        cause_set_id=CONF.get_id("cause_set"),
    )
    scl = pd.read_csv("FILEPATH")
    modeled_causes = scl.loc[scl[int_cause] == 1, "cause_id"].unique().tolist()
    parent_level = chh.loc[chh["cause_id"] == parent_cause_id, "level"].iloc[0]
    if parent_level > 0:
        all_child_causes = get_all_child_cause_ids(parent_cause_id, chh)
    else:
        all_child_causes = chh.loc[
            (chh["parent_id"] == 294) & (chh["level"] == 1) & (chh["secret_cause"] == 0), "cause_id"
        ].unique()

    child_causes = chh.loc[
        (chh["cause_id"].isin(all_child_causes))
        & (chh["level"] == parent_level + 1)
        & (chh["secret_cause"] == 0),
        "cause_id",
    ].unique()
    if parent_level == 3:
        child_causes = list(set(child_causes).intersection(set(modeled_causes)))
    elif parent_level <= 2:
        if int_cause == "sepsis":
            model_folder = "01_sepsis"
        else:
            model_folder = "02_infectious_syndrome"
        prev_level = pd.read_parquet("FILEPATH")
        prev_level_agg_causes = prev_level["cause_id"].unique().tolist()
        if parent_level == 2:
            modeled_l3_causes = (
                scl.loc[(scl["cause_id"].isin(modeled_causes)) & (scl["level"] == 3), "cause_id"]
                .unique()
                .tolist()
            )
            prev_level_agg_causes += modeled_l3_causes
        child_causes = list(set(child_causes).intersection(set(prev_level_agg_causes)))

    cause_dfs = []
    for cause_id in child_causes:
        logger.debug(f"Working on child cause {cause_id}")
        cause_df = (
            McauseResult(
                int_cause=int_cause,
                end_product=end_product,
                process="age_loc_aggregation",
                year_id=year_id,
                cause_id=cause_id,
                description=description,
                conf=CONF,
            )
            .read_results()
            .query("age_group_id != 27")
        )
        cause_dfs.append(cause_df)
    df = pd.concat(cause_dfs, ignore_index=True, sort=True)
    df["cause_id"] = parent_cause_id
    df = df.groupby(DEM_COLS, as_index=False)[EST_COLS].sum()
    if end_product in ["mortality", "incidence"]:
        logger.debug("Creating age standardized rates")
        age_std_df = create_age_standardized_rates(df, EST_COLS, DEM_COLS)
        df = pd.concat([age_std_df, df], ignore_index=True, sort=True)
        df[DEM_COLS] = df[DEM_COLS].astype(int)
    report_duplicates(df, DEM_COLS)
    logger.debug("Saving output")
    McauseResult(
        int_cause=int_cause,
        end_product=end_product,
        process="cause_aggregation",
        year_id=year_id,
        cause_id=parent_cause_id,
        description=description,
        conf=CONF,
    ).write_results(df)


if __name__ == "__main__":
    description = str(sys.argv[1])
    year_id = str(sys.argv[2])
    int_cause = str(sys.argv[3])
    end_product = str(sys.argv[4])
    parent_cause_id = int(sys.argv[5])
    main_cause_agg(description, year_id, int_cause, end_product, parent_cause_id)
