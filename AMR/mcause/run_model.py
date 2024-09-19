import argparse
import getpass
import logging
import math
import multiprocessing
import os
import random
import re
import shlex
import subprocess
import sys
import time
from multiprocessing import Pool, current_process

import pandas as pd
from tqdm import tqdm

from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import (
    add_cause_metadata,
    add_location_metadata,
    create_age_bins,
    get_all_child_cause_ids,
    get_cod_ages,
    get_current_cause_hierarchy,
    prep_child_to_available_parent_map,
    pretty_print,
)
from cod_prep.utils import print_log_message, report_if_merge_fail
from mcod_prep.utils.causes import get_int_cause_hierarchy
from mcod_prep.utils.covariates import merge_covariate
from mcod_prep.utils.mcause_io import McauseResult, get_mcause_data, makedirs_safely, setup_logger
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from mcod_prep.utils.nids import get_datasets

CONF = Configurator("standard")
NUM_WORKER = 26


def fix_sepsis(df, int_cause):
    assert int_cause in ["explicit_sepsis", "implicit_sepsis", "sepsis"]
    sepsis_type = int_cause.split("_")[0]
    if sepsis_type == "sepsis":
        df[int_cause] = (df["sepsis"].isin(["explicit", "implicit"])) * 1
    else:
        df[int_cause] = (df["sepsis"] == sepsis_type) * 1


def get_all_scl_causes(int_cause):
    chh = get_current_cause_hierarchy(
        cause_set_version_id=CONF.get_id("cause_set_version"),
        cause_set_id=CONF.get_id("cause_set"),
    )

    syn_cause_map = pd.read_csv("FILEPATH")
    parents = syn_cause_map.loc[syn_cause_map[int_cause] == 1, "cause_id"].unique().tolist()

    full_cause_list = get_all_child_cause_ids(parents, chh)
    full_cause_list = list(set(full_cause_list))

    missing_elements = [x for x in parents if x not in full_cause_list]
    missing_elements = list(set(missing_elements))
    print_log_message(
        f"The full cause list is missing cause_ids from the syndrome cause list: {missing_elements}"
    )
    full_cause_list += missing_elements

    child_parent_dict = prep_child_to_available_parent_map(
        chh, parents, as_dict=True, expect_all=False
    )

    return full_cause_list, child_parent_dict, parents


def agg_child_causes(df, int_cause):
    deaths_before = df.deaths.sum()
    full_cause_list, child_parent_dict, parents = get_all_scl_causes(int_cause)
    df.loc[
        (df["level_2"].isin(full_cause_list)) & (~df["level_2"].isin(parents)), "remap_child_cause"
    ] = 1
    df["remap_child_cause"].fillna(0, inplace=True)

    temp_remap_df = df.loc[df["remap_child_cause"] == 1]
    temp_nochange_df = df.loc[df["remap_child_cause"] == 0]
    temp_remap_df["level_2"] = temp_remap_df["level_2"].apply(lambda x: child_parent_dict.get(x, x))

    temp_remap_df.drop(columns=["cause_id", "level_1"], inplace=True)
    print_log_message(
        "Pulling in cause hierarchy to reassign 'level_1' and 'cause_id' to remapped children."
    )
    cause_meta_df = get_int_cause_hierarchy(
        int_cause,
        force_rerun=False,
        block_rerun=True,
        cache_dir="standard",
        cause_set_version_id=CONF.get_id("cause_set_version"),
    )[["cause_id", "level_1", "level_2"]]

    temp_remap_df = temp_remap_df.merge(cause_meta_df, how="left", on="level_2")
    df = pd.concat([temp_remap_df, temp_nochange_df], ignore_index=True)

    drop_cols = ["remap_child_cause", "yll_only", "yld_only"]
    df.drop(columns=[x for x in drop_cols if x in df.columns], inplace=True)
    group_cols = [x for x in df.columns if x not in ["deaths"]]
    print(f"RIGHT BEFORE GROUPBY {df.deaths.sum()}")
    df_grouped = df.groupby(group_cols, as_index=False)["deaths"].sum()
    deaths_after = df_grouped.deaths.sum()

    if deaths_before != deaths_after:
        print_log_message(
            f"FLAG! Total number of deaths changed during child remapping from {deaths_before} to {deaths_after}."
        )
    return df_grouped


def collapse_across_value_cols(df, end_product, int_cause, id_cols):
    value_cols = ["successes", "failures"]
    denominator = "cases"
    if end_product == "attributable_burden":
        df["successes"] = (df[int_cause] != 1) * df["deaths"]
        df["cases"] = (df[int_cause] != 1) * df["admissions"]
    else:
        df["successes"] = (df[int_cause] == 1) * df["deaths"]
        if end_product == "incidence":
            df["cases"] = (df[int_cause] == 1) * df["admissions"]
        elif end_product == "mortality":
            denominator = "deaths"
    if denominator == "cases":
        df = df[df["cases"] > 0]
    df["failures"] = df[denominator] - df["successes"]
    value_cols += [denominator]
    print_log_message(f"Collapsing {value_cols} across {id_cols}")
    df = df.groupby(id_cols, as_index=False)[value_cols].sum()
    df["obs_fraction"] = df["successes"] / df[denominator]
    return df


def format_for_model(
    df,
    int_cause,
    end_product,
    age_group_ids,
    id_cols=["year_id", "sex_id", "age_group_id", "location_id", "level_1", "level_2"],
):

    if int_cause in ["L2_genital_infection"]:
        keep_ages = get_cod_ages()["age_group_id"].unique().tolist()
        keep_ages = [x for x in keep_ages if x > 6]
        dropna = True
    else:
        keep_ages = get_cod_ages()["age_group_id"].unique().tolist()
        dropna = False

    drop_rows = (
        (df["sex_id"] == 9)
        | ~(df["age_group_id"].isin(keep_ages))
        | (df["cause_id"].isin([919, 744, 743]))
    )
    df = df[~drop_rows]

    df = create_age_bins(df, age_group_ids, dropna=dropna)

    df = add_cause_metadata(
        df,
        ["yld_only", "yll_only"],
        block_rerun=False,
        force_rerun=True,
        cause_set_version_id=CONF.get_id("cause_set_version"),
    )
    if end_product == "mortality":
        df = df[(df["deaths"] > 0) & (df["yld_only"] != 1)]
        if int_cause in CONF.get_id("inf_syns") + ["implicit_sepsis", "explicit_sepsis", "sepsis"]:
            df = agg_child_causes(df, int_cause)
    else:
        df = df[(df["admissions"] > 0) & (df["yll_only"] != 1)]
    if int_cause in ["implicit_sepsis", "explicit_sepsis", "sepsis"]:
        fix_sepsis(df, int_cause)
    assert int_cause in df.columns, "intermediate cause is missing!"
    assert set(df[int_cause].unique()) == {0, 1}, f"expecting {int_cause} column to be 0 or 1"
    df = collapse_across_value_cols(df, end_product, int_cause, id_cols)
    df[id_cols] = df[id_cols].astype(int)
    return df


def merge_nested_cause_levels(df, int_cause):
    print_log_message("Pulling in cause hierarchy")
    cause_meta_df = get_int_cause_hierarchy(
        int_cause,
        force_rerun=False,
        block_rerun=True,
        cache_dir="standard",
        cause_set_version_id=CONF.get_id("cause_set_version"),
    )[["cause_id", "level_1", "level_2"]]

    df = df.merge(cause_meta_df, how="left", on="cause_id")
    return df


def pull_covariates(int_cause):
    covariates_df = pd.read_csv("FILEPATH").query("int_cause == @int_cause")
    assert len(covariates_df) == 1
    return covariates_df["covariates"].str.split(",").iloc[0]


def write_model_input_data(df, int_cause, end_product, description, output_dir, age_group_id=None):
    makedirs_safely(output_dir)
    print_log_message(f"Writing model input file to {output_dir}")
    if age_group_id:
        df = df.query(f"age_group_id == {age_group_id}")
    df.to_csv("FILEPATH", index=False)


def launch_modelworker(args):
    int_cause, end_product, description, oosv, age_group_id = args
    log_path = "FILEPATH"

    output_dir = "FILEPATH"

    diag_dir = "FILEPATH"
    log_dir = "FILEPATH"
    in_sample_dir = "FILEPATH"
    oo_sample_dir = "FILEPATH"
    makedirs_safely(in_sample_dir)
    makedirs_safely(oo_sample_dir)
    makedirs_safely(log_dir)

    logger = setup_logger(
        f"{int_cause}_{end_product}_{description}_{age_group_id}_modelworker", log_dir
    )

    worker = "FILEPATH"
    logger.debug(f"Working on age {age_group_id}")

    cmd = [
        "FILEPATH",
        "-s",
        worker,
        int_cause,
        output_dir,
        diag_dir,
        str(CONF.get_id("release")),
        str(oosv),
    ]

    cmd_str = shlex.join(cmd)

    try:
        with subprocess.Popen(
            cmd_str,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            text=True,
            bufsize=1,
            universal_newlines=True,
        ) as proc:
            for line in proc.stdout:
                logger.info(f"R script output for age {age_group_id}: {line}")
            for line in proc.stderr:
                logger.error(f"R script error for age {age_group_id}: {line}")

            proc.wait()
            if proc.returncode != 0:
                logger.error(
                    f"R script failed with return code {proc.returncode} for age {age_group_id}"
                )

    except Exception as e:
        logger.error(f"Exception occurred during subprocess call for age {age_group_id}:\n{e}")


def get_SOURCE(int_cause):
    df = get_mcause_data(
        phase="format_map",
        sub_dirs=int_cause,
        assert_all_available=True,
        force_rerun=True,
        block_rerun=False,
        is_active=True,
        source="SOURCE",
    )
    group_cols = list(set(df.columns) - set(["code_id", "deaths"]))
    return (
        df.query("age_group_id != 6 & cause_id != 743")
        .groupby(group_cols, as_index=False)["deaths"]
        .sum()
    )


def prep_input_data(description, int_cause, end_product, project_id, age_group_ids):
    data_kwargs = {
        "phase": "redistribution",
        "sub_dirs": int_cause,
        "assert_all_available": False,
        "force_rerun": True,
        "block_rerun": False,
        "is_active": True,
        "data_type_id": [9, 3, 13],
        "year_id": list(range(1980, 2050)),
        "project_id": project_id,
    }
    if end_product in ["incidence", "attributable_burden"]:
        data_kwargs.update({"data_type_id": 3})
    print_log_message("Pulling training data")
    df = get_mcause_data(**data_kwargs)
    if end_product == "mortality":
        NUM_LIST = ["NUMBERS"]
        df = df.query(f"nid not in {NUM_LIST}")
        print_log_message("Appending SOURCE")
        df = pd.concat([df, get_SOURCE(int_cause)], sort=True)
    print_log_message("Merging nested cause levels")
    df = merge_nested_cause_levels(df, int_cause)
    print_log_message("Formatting for model")
    df = format_for_model(df, int_cause, end_product, age_group_ids)
    print_log_message("Pulling covariates")
    covariates = pull_covariates(int_cause)
    for covariate in covariates:
        print_log_message("Merging on {}".format(covariate))
        df = merge_covariate(df, covariate)
    assert df.notnull().values.all()
    if "by_age" in description:
        for age_group_id in age_group_ids:
            print_log_message("Saving Model Input")
            output_dir_by_age = "FILEPATH"
            write_model_input_data(
                df, int_cause, end_product, description, output_dir_by_age, age_group_id
            )
    else:
        output_dir = str("FILEPATH")
        write_model_input_data(df, int_cause, end_product, description, output_dir)


def main(description, int_cause, end_product, project_id, oosv, age_group_ids):
    print_log_message("Prepping input data")
    prep_input_data(description, int_cause, end_product, project_id, age_group_ids)

    if "by_age" in description:
        print_log_message("Launching model by age")
        with multiprocessing.Manager() as manager:
            with Pool(processes=NUM_WORKER) as pool:
                for i in tqdm(
                    pool.imap_unordered(
                        launch_modelworker,
                        [
                            [int_cause, end_product, description, oosv, age_group_id]
                            for age_group_id in age_group_ids
                        ],
                    )
                ):
                    pass
    else:
        print("Parallelizing not by age is not available yet")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format data and launch model.")
    parser.add_argument("description", type=str)
    parser.add_argument("int_cause", type=str)
    parser.add_argument("end_product", type=str)
    parser.add_argument("project_id", type=int)
    parser.add_argument("oosv", type=int)
    parser.add_argument("--age_group_ids", action="append", type=int)
    args = parser.parse_args()
    main(**vars(args))
