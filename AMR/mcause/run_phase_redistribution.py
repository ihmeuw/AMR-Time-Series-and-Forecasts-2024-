import argparse

import numpy as np

from cod_prep.claude.configurator import Configurator
from cod_prep.claude.redistribution import Redistributor
from cod_prep.claude.run_phase_redistribution import cause_map_post_rd
from cod_prep.downloaders import (
    get_cause_map,
    get_cod_ages,
    get_current_cause_hierarchy,
    get_current_location_hierarchy,
    get_package_map,
)
from cod_prep.utils import CodSchema, clean_icd_codes, print_log_message
from mcod_prep.mcod_mapping import MCoDMapper
from mcod_prep.utils.causes import get_all_related_syndromes
from mcod_prep.utils.logging import ymd_timestamp
from mcod_prep.utils.mcause_io import get_phase_output, write_phase_output

CONF = Configurator("standard")
CACHE_DIR = "FILEPATH"


def reassign_syndrome(df, code_meta_df, int_cause, code_system_id, code_map_version_id):
    mapper = MCoDMapper(
        int_cause,
        drop_p2=True,
        code_system_id=code_system_id,
        code_map_version_id=code_map_version_id,
    )
    int_cause_map = mapper.prep_int_cause_map()
    code_meta_df = code_meta_df.assign(
        cause=lambda d: clean_icd_codes(d["value"], remove_decimal=True)
    )
    df = df.merge(code_meta_df.loc[:, ["code_id", "cause"]])
    df = mapper.map_cause_codes(df, int_cause_map, int_cause, cols_to_map=["cause"])
    mapped_col = "cause_" + int_cause
    df[mapped_col] = df[mapped_col].fillna(MCoDMapper.remainder_syndromes[0])
    if int_cause == "infectious_syndrome":
        df.loc[~df[mapped_col].isin(MCoDMapper.remainder_syndromes), int_cause] = df[mapped_col]
    else:

        if "_comm" not in int_cause and "_hosp" not in int_cause:
            target_syndromes = get_all_related_syndromes(int_cause)

            df.loc[df[mapped_col].isin(target_syndromes), int_cause] = 1

            df.loc[
                ~df[mapped_col].isin(target_syndromes)
                & ~df[mapped_col].isin(MCoDMapper.remainder_syndromes),
                int_cause,
            ] = 0
        if "_comm" in int_cause:
            target_syndromes = get_all_related_syndromes(int_cause.replace("_comm", ""))

            df.loc[df[mapped_col].isin(target_syndromes), int_cause] = 1

            df.loc[~df[mapped_col].isin(target_syndromes), int_cause] = 0

        if "_hosp" in int_cause:
            target_syndromes = get_all_related_syndromes(int_cause.replace("_hosp", ""))

            df.loc[df[mapped_col].isin(target_syndromes), int_cause] = 0

            df.loc[
                ~df[mapped_col].isin(MCoDMapper.remainder_syndromes)
                & ~df[mapped_col].isin(target_syndromes),
                int_cause,
            ] = 0

    df = df.drop([mapped_col, "cause"], axis="columns")

    return df


def trim_icd_len(row):
    if row["code_system_id"] == 6 and row["value"][0] in ["E", "V"] and len(row["value"]) > 5:
        return row["value"][:5]
    elif len(row["value"]) > 4:
        return row["value"][:4]
    else:
        return row["value"]


def prep_trim_dict(cause_map, missing_codeids):

    df = cause_map.copy()
    df = df[["code_system_id", "value", "code_id"]]
    df["value"] = clean_icd_codes(df["value"], remove_decimal=True)

    icd10_mask = (df["code_system_id"] == 1) & df["value"].str.match("^[A-Z][0-9A-Z]*$")
    icd9_mask = (df["code_system_id"] == 6) & df["value"].str.match("^[EV0-9][0-9A-Z]*$")
    df = df[icd10_mask | icd9_mask]

    df["trimmed_value"] = df.apply(trim_icd_len, axis=1)

    df["too_long"] = (
        (
            (df["code_system_id"] == 6)
            & df["value"].str.startswith(("E", "V"))
            & (df["value"].str.len() > 5)
        )
        | (
            (df["code_system_id"] == 6)
            & ~df["value"].str.startswith(("E", "V"))
            & (df["value"].str.len() > 4)
        )
        | ((df["code_system_id"] == 1) & (df["value"].str.len() > 4))
    )
    code_id_map = df.loc[~df["too_long"]].set_index("value")["code_id"].to_dict()
    df["new_code_id"] = df["trimmed_value"].map(code_id_map)

    df.loc[
        df["trimmed_value"].str.contains("E81") & df["new_code_id"].isna(), "new_code_id"
    ] = 157003.0
    df.loc[
        df["trimmed_value"].str.contains("E82") & df["new_code_id"].isna(), "new_code_id"
    ] = 103640.0
    df.loc[df["value"].str.contains("T89") & df["new_code_id"].isna(), "new_code_id"] = 38783.0

    df.loc[df["value"].str.contains("2956") & df["new_code_id"].isna(), "new_code_id"] = 5856.0

    nan_values = df.loc[df["new_code_id"].isna(), "value"]

    if nan_values.any():
        raise ValueError(
            f"NaNs in 'new_code_id' for the following icd codes: {nan_values.tolist()}"
        )

    trim_dict = df.set_index("code_id")["new_code_id"].to_dict()

    trim_dict = {k: v for k, v in trim_dict.items() if k in missing_codeids}

    return trim_dict


def run_phase(
    df,
    csvid,
    lsvid,
    cmvid,
    csid,
    remove_decimal,
    value_col,
    int_cause,
):
    read_file_cache_options = {
        "block_rerun": True,
        "cache_dir": CACHE_DIR,
        "force_rerun": False,
        "cache_results": False,
    }
    cause_map = get_cause_map(code_map_version_id=cmvid, **read_file_cache_options)
    pkg_map = get_package_map(csid)
    lhh = get_current_location_hierarchy(location_set_version_id=lsvid, **read_file_cache_options)
    age_meta_df = get_cod_ages(**read_file_cache_options)
    cause_meta_df = get_current_cause_hierarchy(csvid, **read_file_cache_options)
    col_meta = {
        c: {"col_type": "demographic"}
        for c in [
            "died",
            "sepsis",
            "code_system_id",
            "pathogen_from_cause",
            "cause_cross_tabulate",
            int_cause,
        ]
    }

    col_meta["admissions"] = {"col_type": "value"}
    orig_sum = df[value_col].sum()

    pkg_map["code_id"], df["code_id"] = pkg_map["code_id"].astype(int), df["code_id"].astype(int)

    missing_codeids = set(df.loc[df["cause_id"] == 743, "code_id"]) - set(pkg_map["code_id"])

    if len(missing_codeids) > 0:
        trim_dict = prep_trim_dict(cause_map, missing_codeids)
        df.loc[df["code_id"].isin(missing_codeids), "code_id"] = df["code_id"].map(trim_dict)

        missing = missing_codeids - set(trim_dict.keys())
        if len(missing) > 0:
            raise ValueError(
                f"There are code_ids in the data that were not in the trim_dict: {missing.tolist()}"
            )

    df.drop(columns="cause_id", inplace=True)

    proportional_cols = {{"location_id": "admin1_or_above_id"}.get(c, c) for c in df.columns} & (
        set(CONF.get_id("potential_rd_proportional_cols")) | {int_cause}
    )

    df = Redistributor(
        conf=CONF,
        code_system_id=csid,
        col_meta=col_meta,
        proportional_cols=list(proportional_cols),
        loc_meta_df=lhh,
        age_meta_df=age_meta_df,
        cause_meta_df=cause_meta_df,
        cause_map=cause_map,
        remove_decimal=remove_decimal,
    ).run_redistribution(df)
    if int_cause in MCoDMapper.infectious_syndromes:
        df = reassign_syndrome(
            df,
            code_meta_df=cause_map,
            int_cause=int_cause,
            code_system_id=csid,
            code_map_version_id=cmvid,
        )

    df = cause_map_post_rd(df, cause_map, col_meta=col_meta)

    post_sum = df[value_col].sum()
    before_after_text = (
        f"Before GC redistribution: {orig_sum}\n" f"After GC redistribution: {post_sum}"
    )
    if not np.isclose(orig_sum, post_sum):
        raise AssertionError(f"Deaths not close.\n{before_after_text}")
    else:
        print_log_message(before_after_text)
    return df


def main(
    nid,
    extract_type_id,
    csvid,
    lsvid,
    csid,
    cmvid,
    remove_decimal,
    data_type_id,
    int_cause,
):

    raw_df = get_phase_output("format_map", nid, extract_type_id, sub_dirs=int_cause)

    fatal_df = raw_df.loc[raw_df["deaths"] > 0]

    if len(fatal_df) > 0:

        if "admissions" in fatal_df.columns:
            fatal_df.drop(columns=["admissions"], inplace=True)

        value_col = "deaths"

        if fatal_df.eval("cause_id == 743").any():
            print_log_message("Running redistribution")

            if "cause_cross_tabulate" in fatal_df.columns:
                fatal_df.drop(columns=["cause_cross_tabulate"], inplace=True)

            fatal_df = run_phase(
                fatal_df, csvid, lsvid, cmvid, csid, remove_decimal, value_col, int_cause
            )

            if int_cause == "cross_tabulate":
                fatal_df["cause_cross_tabulate"] = fatal_df["cause_id"]
        else:
            print_log_message("No redistribution to do.")

        df = fatal_df
        group_cols = list(set(df.columns) - set(["code_id", "deaths"]))
        df = df.groupby(group_cols, as_index=False).agg({value_col: "sum"})

    else:
        print_log_message(f"No deaths in source {nid} - no redistribution to do.")
        return

    write_phase_output(
        df, "redistribution", nid, extract_type_id, ymd_timestamp(), sub_dirs=int_cause
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parameters for running redistribution")
    parser.add_argument("nid", type=int)
    parser.add_argument("extract_type_id", type=int)
    parser.add_argument("csvid", type=int)
    parser.add_argument("lsvid", type=int)
    parser.add_argument("csid", type=int)
    parser.add_argument("cmvid", type=int)
    parser.add_argument("remove_decimal", type=bool)
    parser.add_argument("data_type_id", type=int)
    parser.add_argument("int_cause", type=str)
    args = parser.parse_args()

    main(**vars(args))
