import sys
from builtins import str
from importlib import import_module

import pandas as pd

from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import get_cod_ages, get_current_location_hierarchy
from cod_prep.downloaders.causes import get_all_related_causes, get_current_cause_hierarchy
from cod_prep.utils import print_log_message, report_duplicates
from cod_prep.utils.schema import ColType
from mcod_prep.mcod_mapping import MCoDMapper
from mcod_prep.utils.age_sex_split import AgeSexSplitter
from mcod_prep.utils.logging import print_memory_timestamp, ymd_timestamp
from mcod_prep.utils.mcause_io import makedirs_safely, write_phase_output
from mcod_prep.utils.unmodeled_demographics_remap import AgeSexRemap

CONF = Configurator()
PROCESS_DIR = "FILEPATH"
DIAG_DIR = "FILEPATH"
ID_COLS = [
    "year_id",
    "sex_id",
    "age_group_id",
    "location_id",
    "cause_id",
    "code_id",
    "nid",
    "extract_type_id",
]


def get_drop_part2(int_cause, source):
    if int_cause in MCoDMapper.infectious_syndromes:
        return True
    elif int_cause == "cross_tabulate":
        return True
    elif int_cause == "custom":
        return False
    else:
        df = pd.read_csv("FILEPATH", dtype={"drop_p2": bool})
        df = df[df["int_cause"] == int_cause]
        report_duplicates(df, "int_cause")
        try:
            return df["drop_p2"].iloc[0]
        except IndexError:
            print(f"{int_cause} was not found in {CONF.get_resource('drop_p2')}")


def get_formatting_method(source, data_type_id, year, drop_p2, nid):
    clean_source = "clean_" + source.lower()
    if source != "SOURCE":
        args = [year, drop_p2]
    else:
        args = [year, drop_p2, nid]
    try:
        formatting_method = getattr(import_module("FILEPATH", f"{clean_source}"))
    except AttributeError:
        print(f"No formatting method found! Check module & main function are named clean_{source}")
    return formatting_method, args


def drop_non_mcause(df, int_cause):
    chain_cols = [x for x in df.columns if ("multiple_cause_" in x) and ("pII" not in x)]
    cause_cols = ["cause"] + chain_cols
    unique_col_num = 2
    assert (
        len(cause_cols) >= unique_col_num
    ), f"there need to be at least {unique_col_num} cause columns"

    df["unique_codes"] = df[cause_cols].apply(lambda x: set([x for x in x if x != "0000"]), axis=1)
    df["unique_code_len"] = df["unique_codes"].apply(lambda x: len(x))

    row_n = len(df)

    df_dedup = df.loc[df["unique_code_len"] >= unique_col_num]
    dedup_row_n = len(df_dedup)
    perc_dropped = ((row_n - dedup_row_n) / row_n) * 100
    print(f"Dropped {perc_dropped}% rows due to not having 2+ unique codes")

    df_dedup.drop(["unique_codes", "unique_code_len"], axis=1, inplace=True)

    return df_dedup


def collect_cause_specific_diagnostics(df, acause_list):
    cause_meta_df = get_current_cause_hierarchy(force_rerun=True, block_rerun=False)
    if type(acause_list) != list:
        acause_list = [acause_list]
    df_list = []
    for acause in acause_list:
        diag_df = df.loc[df["cause_id"].isin(get_all_related_causes(acause, cause_meta_df))]
        df_list.append(diag_df)
    diag_df = pd.concat(df_list, ignore_index=True)
    return diag_df


def injuries_diagnostic(df):
    patterns = df["pattern"]
    df = df[(patterns.isin(patterns[patterns.duplicated()])) & (df.pattern != "")]
    return df


def get_id_value_cols(df, int_cause, data_type_id, inj_diag, drop_p2):
    group_cols = ID_COLS + [int_cause]
    if int_cause in ["x59", "y34"]:
        group_cols += ["pattern"]
        if inj_diag:
            group_cols += ["pII_ncodes", "pII_in_ncodes"]
    elif int_cause == "infectious_syndrome":
        group_cols += [
            "sepsis",
            "pathogen_from_cause",
            "cause_infectious_syndrome",
            "all_possible_syndromes",
            "most_severe",
        ]
    elif int_cause == "cross_tabulate":
        group_cols += ["cause_cross_tabulate"]
    elif int_cause == "standard_map":
        group_cols += ["code_system_id"]
        group_cols += [x for x in df.columns if "standard_map" in x]
    elif int_cause == "flu_type":
        group_cols += [
            "cause_flu_type",
            "infectious_syndrome",
            "cause_infectious_syndrome",
        ]
    elif int_cause == "pathogen":
        group_cols = [x for x in df.columns if x not in ["admissions", "deaths"]]
    if not drop_p2:
        group_cols += ["pII_" + int_cause]
    group_cols = list(set(group_cols))
    if data_type_id in [3, 17]:
        value_cols = ["admissions", "deaths"]
    else:
        value_cols = ["deaths"]
    return group_cols, value_cols


def run_pipeline(
    year,
    source,
    int_cause,
    code_system_id,
    code_map_version_id,
    nid,
    extract_type_id,
    data_type_id,
    inj_diag=False,
):
    drop_p2 = get_drop_part2(int_cause, source)

    print_log_message("Prepping data")
    formatting_method, args = get_formatting_method(source, data_type_id, year, drop_p2, nid)
    df = formatting_method(*args)

    print_log_message("Dropping rows without multiple cause")
    df = drop_non_mcause(df, int_cause)
    assert len(df) > 0, "No multiple cause data here!"

    print_log_message("Mapping data")
    Mapper = MCoDMapper(
        int_cause, drop_p2, code_system_id=code_system_id, code_map_version_id=code_map_version_id
    )
    df = Mapper.get_computed_dataframe(df)

    if (int_cause in MCoDMapper.infectious_syndromes) and (int_cause != "pathogen"):
        print_log_message("Mapping data again to get sepsis-related observations")
        Mapper = MCoDMapper(
            "sepsis",
            drop_p2,
            code_system_id=code_system_id,
            code_map_version_id=code_map_version_id,
        )
        df = Mapper.get_computed_dataframe(df, map_underlying_cause=False)

        if int_cause not in ["infectious_syndrome"]:
            print_log_message("Dropping non-sepsis related observations")
            df = df.loc[
                (df["sepsis"] != "no_sepsis")
                & (df["infectious_syndrome"] != "unspecified_and_none")
            ]
            assert len(df) > 0, "No sepsis-related rows!"
    elif int_cause == "flu_type":
        print_log_message("Mapping data again to get infectious syndromes")
        Mapper = MCoDMapper(
            "infectious_syndrome",
            drop_p2,
            code_system_id=code_system_id,
            code_map_version_id=code_map_version_id,
        )
        df = Mapper.get_computed_dataframe(df, map_underlying_cause=False)

    if inj_diag and (int_cause in ["x59", "y34"]):
        ddf = injuries_diagnostic(df)
        outdir = DIAG_DIR.format(nid=nid, extract_type_id=extract_type_id, int_cause=int_cause)
        makedirs_safely(outdir)
        ddf.to_csv("FILEPATH", index=False)

    group_cols, value_cols = get_id_value_cols(df, int_cause, data_type_id, inj_diag, drop_p2)

    if int_cause != "standard_map":
        print_memory_timestamp(df, f"Collapsing {value_cols} across {group_cols}")
        df = df.groupby(group_cols, as_index=False)[value_cols].sum()
    else:
        print_memory_timestamp(df, "Standard map does not collapse cols")
        df.drop(columns=["standard_map"], inplace=True)

    print_memory_timestamp(df, "Filtering cause-age-sex restrictions")
    cause_set_version_id = CONF.get_id("cause_set_version")
    column_metadata = {
        **{group_col: {"col_type": ColType.DEMOGRAPHIC} for group_col in group_cols},
        **{value_col: {"col_type": ColType.VALUE} for value_col in value_cols},
    }

    Remapper = AgeSexRemap(
        code_system_id,
        cause_set_version_id,
        int_cause,
        collect_diagnostics=False,
        verbose=True,
        column_metadata=column_metadata,
    )
    if int_cause != "standard_map":
        df = Remapper.get_computed_dataframe(df)
    else:
        df = Remapper.get_computed_dataframe(df, int_cause)

    if int_cause in ["flu_type", "infectious_syndrome"]:
        group_cols = list(set(group_cols) - {"code_id"})
        df = df.groupby(group_cols, as_index=False)[value_cols].sum()

    Splitter = AgeSexSplitter(verbose=True, collect_diagnostics=True)

    df_age_ids = df["age_group_id"].unique().tolist()
    age_meta_df = get_cod_ages()
    cod_age_ids = age_meta_df["age_group_id"].unique().tolist()

    df_sex_ids = df["sex_id"].unique().tolist()
    if any(age not in cod_age_ids for age in df_age_ids) or any(
        sex in [3, 9] for sex in df_sex_ids
    ):
        print_log_message("Performing age-sex splitting")
        lhh = get_current_location_hierarchy()

        age_sex_df = Splitter.get_computed_dataframe(df, lhh)
    else:
        print_log_message("No age-sex splitting needed")

        age_sex_df = df

    return df, age_sex_df


def main(
    year, source, int_cause, code_system_id, code_map_version_id, nid, extract_type_id, data_type_id
):
    df, age_sex_df = run_pipeline(
        year,
        source,
        int_cause,
        code_system_id,
        code_map_version_id,
        nid,
        extract_type_id,
        data_type_id,
    )
    print_log_message(f"Writing nid {nid}, extract_type_id {extract_type_id}")
    if int_cause == "pathogen":
        write_phase_output(
            df,
            "format_map_raw",
            nid,
            extract_type_id,
            ymd_timestamp(),
            sub_dirs=int_cause,
            file_format="parquet",
        )

        write_phase_output(
            age_sex_df,
            "format_map",
            nid,
            extract_type_id,
            ymd_timestamp(),
            sub_dirs=int_cause,
            file_format="parquet",
        )
    else:
        write_phase_output(
            df, "format_map_raw", nid, extract_type_id, ymd_timestamp(), sub_dirs=int_cause
        )

        write_phase_output(
            age_sex_df, "format_map", nid, extract_type_id, ymd_timestamp(), sub_dirs=int_cause
        )


if __name__ == "__main__":
    year = int(sys.argv[1])
    source = str(sys.argv[2])
    int_cause = str(sys.argv[3])
    code_system_id = int(sys.argv[4])
    code_map_version_id = int(sys.argv[5])
    nid = int(sys.argv[6])
    extract_type_id = int(sys.argv[7])
    data_type_id = int(sys.argv[8])

    main(
        year,
        source,
        int_cause,
        code_system_id,
        code_map_version_id,
        nid,
        extract_type_id,
        data_type_id,
    )
