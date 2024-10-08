import getpass
import logging
import multiprocessing
import os
import re
import time
import warnings
from datetime import datetime, timedelta
from multiprocessing import Pool, current_process
from pathlib import Path, PosixPath
from typing import List

import numpy as np
import pandas as pd

from db_queries import get_ids
from db_tools import ezfuncs

from cod_prep.claude.configurator import Configurator
from cod_prep.utils import cod_timestamp, print_log_message, wrap
from cod_prep.utils.misc import print_log_message, report_duplicates, umask
from cod_prep.utils.read_write_df import read_df, write_df
from mcod_prep.mcod_mapping import MCoDMapper
from mcod_prep.utils.nids import get_datasets

CONF = Configurator("standard")


def makedirs_safely(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError as e:
        if e.errno == errno.EEXIST:
            print(
                "Process could not create directory {} because it already"
                "existed, probably due to race condition, no error,"
                "continuing".format(directory)
            )
            pass
        else:
            print("Process could not create directory {}, " "re-raising".format(directory))
            raise


def get_overlapping_nid_extracts(df, table, conn_def="", primary_keys=["nid", "extract_type_id"]):
    if "project_id" in primary_keys:
        db_nid_extracts = ezfuncs.query(
            """
            SELECT project_id, nid, extract_type_id
            FROM cod.{tbl}
        """.format(
                tbl=table
            ),
            conn_def=conn_def,
        )
    else:
        db_nid_extracts = ezfuncs.query(
            """
            SELECT nid, extract_type_id
            FROM cod.{tbl}
        """.format(
                tbl=table
            ),
            conn_def=conn_def,
        )
    db_nid_extracts = db_nid_extracts.drop_duplicates()
    nid_extracts_df = df[primary_keys].drop_duplicates()
    overlap = db_nid_extracts.merge(nid_extracts_df)
    return overlap


def delete_nid_extracts_from_table(table, nid_extracts, conn_def="prodcod"):
    report_duplicates(nid_extracts, ["nid", "extract_type_id"])
    nid_extract_list = nid_extracts[["nid", "extract_type_id"]].to_records(index=False)

    delete_query = """
    DELETE
    FROM cod.{tbl}
    WHERE (
        {nid_extract_clause}
    )
    """

    nid_extract_clause = ""
    first_nid_ex = "(nid={nid} AND extract_type_id={extract_type_id})"
    next_nid_ex_template = "\n    OR (nid={nid} AND extract_type_id={extract_type_id})"

    first_nid_ex_rec = nid_extract_list[0]
    nid = int(first_nid_ex_rec[0])
    extract_type_id = int(first_nid_ex_rec[1])
    first_nid_ex = first_nid_ex.format(nid=nid, extract_type_id=extract_type_id)
    nid_extract_clause = nid_extract_clause + first_nid_ex
    for nid, extract_type_id in nid_extract_list[1:]:
        nid = int(nid)
        next_nid_ex = next_nid_ex_template.format(nid=nid, extract_type_id=extract_type_id)
        nid_extract_clause = nid_extract_clause + next_nid_ex
    delete_query = delete_query.format(tbl=table, nid_extract_clause=nid_extract_clause)

    engine = ezfuncs.get_engine(conn_def=conn_def)
    conn = engine.connect()
    res = conn.execute(delete_query)
    conn.close()
    print("Successfully deleted {} rows from {}".format(res.rowcount, table))


def delete_project_nid_extracts_from_table(table, proj_nid_extracts, conn_def="prodcod"):
    report_duplicates(proj_nid_extracts, ["project_id", "nid", "extract_type_id"])
    proj_nid_extract_list = proj_nid_extracts[["project_id", "nid", "extract_type_id"]].to_records(
        index=False
    )

    delete_query = """
    DELETE
    FROM cod.{tbl}
    WHERE (
        {proj_nid_extract_clause}
    )
    """

    proj_nid_extract_clause = ""
    first_proj_nid_ex = (
        "(project_id={project_id} AND nid={nid} AND extract_type_id={extract_type_id})"
    )
    next_proj_nid_ex_template = (
        "\n    OR (project_id={project_id} AND nid={nid} AND extract_type_id={extract_type_id})"
    )

    first_proj_nid_ex_rec = proj_nid_extract_list[0]
    project_id = int(first_proj_nid_ex_rec[0])
    nid = int(first_proj_nid_ex_rec[1])
    extract_type_id = int(first_proj_nid_ex_rec[2])
    first_proj_nid_ex = first_proj_nid_ex.format(
        project_id=project_id, nid=nid, extract_type_id=extract_type_id
    )
    proj_nid_extract_clause = proj_nid_extract_clause + first_proj_nid_ex
    for project_id, nid, extract_type_id in proj_nid_extract_list[1:]:
        project_id = int(project_id)
        nid = int(nid)
        next_proj_nid_ex = next_proj_nid_ex_template.format(
            project_id=project_id, nid=nid, extract_type_id=extract_type_id
        )
        proj_nid_extract_clause = proj_nid_extract_clause + next_proj_nid_ex
    delete_query = delete_query.format(tbl=table, proj_nid_extract_clause=proj_nid_extract_clause)

    engine = ezfuncs.get_engine(conn_def=conn_def)
    conn = engine.connect()
    res = conn.execute(delete_query)
    conn.close()
    print("Successfully deleted {} rows from {}".format(res.rowcount, table))


def construct_phase_path(
    phase,
    nid,
    extract_type_id,
    file_format="standard",
    sub_dirs=None,
    launch_set_id=None,
    makedirs=True,
    refresh_id=None,
    project_id=None,
    verbose=False,
):
    assert CONF.config_type in ["mcod", "amr", "mcod_rdp", "mcod_for_covid"]

    if sub_dirs is not None:
        if sub_dirs in CONF.get_id("inf_syns"):
            assert CONF.config_type == "amr"
        elif sub_dirs == "sepsis":
            assert CONF.config_type in ["amr", "mcod", "mcod_rdp"]

    if project_id == None:
        project_id = CONF.get_id("project")

    if verbose == True:
        print(f"Configuring path to the result for {CONF.config_type}, project_id {project_id}")

    avail_formats = ["standard", "search", "h5", "csv", "parquet"]
    assert file_format in avail_formats, "Unsupported file format: {}. " "Choose from {}".format(
        file_format, avail_formats
    )
    file_format_standard = CONF.get_id("file_format")
    if file_format == "standard":
        file_format = file_format_standard

    data_dir = "FILEPATH"
    if refresh_id is not None:
        data_dir = "FILEPATH"
    nid_extract_type_id_dir = "FILEPATH"
    if sub_dirs is not None:
        nid_extract_type_id_dir = "FILEPATH"
    arch_dir = "FILEPATH"

    if makedirs:
        if sub_dirs != "_diagnostic":
            makedirs_safely(arch_dir)
        if sub_dirs is not None:
            makedirs_safely(nid_extract_type_id_dir)

    if launch_set_id is None:
        file_name = "FILEPATH"
    else:
        file_name = "FILEPATH"

    if file_format == "search":
        csv_exists = os.path.exists("FILEPATH")
        h5_exists = os.path.exists("FILEPATH")
        parquet_exists = os.path.exists("FILEPATH")
        file_format = ""
        if csv_exists & h5_exists:
            file_format = file_format_standard
        elif csv_exists:
            file_format = "csv"
        elif h5_exists:
            file_format = "h5"
        elif parquet_exists:
            file_format = "parquet"
        else:
            file_format = file_format_standard
        assert file_format != "", "impossible bug"

    file_name = file_name + "." + file_format

    return file_name


def write_phase_output(
    df,
    phase,
    nid,
    extract_type_id,
    launch_set_id,
    queryable_columns=["cause_id"],
    sub_dirs=None,
    archive=True,
    file_format="standard",
):
    with umask(0o002):
        nid = int(nid)
        extract_type_id = int(extract_type_id)
        outpath = construct_phase_path(
            phase,
            nid,
            extract_type_id,
            file_format=file_format,
            sub_dirs=sub_dirs,
            makedirs=True,
        )

        if archive:
            archpath = construct_phase_path(
                phase,
                nid,
                extract_type_id,
                file_format=file_format,
                sub_dirs=sub_dirs,
                launch_set_id=launch_set_id,
                makedirs=True,
            )
            if file_format == "parquet":
                cause_cols = [
                    x for x in df.columns if "cause" in x and x not in ["cause_id", "code_id"]
                ]
                df[cause_cols] = df[cause_cols].astype(str)
            write_df(df, archpath, queryable_columns)

        if launch_set_id != 0:
            if file_format == "parquet":
                cause_cols = [
                    x for x in df.columns if "cause" in x and x not in ["cause_id", "code_id"]
                ]
                df[cause_cols] = df[cause_cols].astype(str)
            write_df(df, outpath, queryable_columns)


def get_phase_output(
    phase,
    nid,
    extract_type_id,
    launch_set_id=None,
    where_filter=None,
    exec_function=None,
    exec_function_args=[],
    sub_dirs=None,
    refresh_id=None,
    usecols=None,
    project_id=None,
):
    filepath = construct_phase_path(
        phase,
        nid,
        extract_type_id,
        launch_set_id=launch_set_id,
        sub_dirs=sub_dirs,
        file_format="search",
        makedirs=False,
        refresh_id=refresh_id,
        project_id=project_id,
    )
    max_tries = 3
    num_tries = 0
    df = None
    while df is None:
        try:
            df = read_df(filepath, where_filter=where_filter, usecols=usecols)
        except IOError:
            if num_tries < max_tries:
                num_tries += 1
                time.sleep(2)
                continue
            else:
                raise
    if exec_function is not None:
        df = exec_function(df, *exec_function_args)
    return df


def check_output_exists(
    phase,
    nid,
    extract_type_id,
    launch_set_id=None,
    sub_dirs=None,
    refresh_id=None,
    project_id=None,
):
    filepath = construct_phase_path(
        phase,
        nid,
        extract_type_id,
        sub_dirs=sub_dirs,
        launch_set_id=launch_set_id,
        makedirs=False,
        refresh_id=refresh_id,
        project_id=project_id,
    )
    filepath = filepath[:-4]
    hdf = filepath + ".h5"
    csv = filepath + ".csv"
    parquet = filepath + ".parquet"
    filepaths = [csv, hdf, parquet]
    exists = []
    for filepath in filepaths:
        if os.path.exists(filepath):
            size = os.path.getsize(filepath)
            if size > 0:
                exists.append(True)
            else:
                exists.append(False)
    if sum(exists) > 0:
        return True
    else:
        return False


def write_to_db_nid_table(
    df,
    table,
    conn_def="",
    replace=False,
    primary_keys=["nid", "extract_type_id"],
):
    df = df.copy()

    assert "mcause" in table

    if "project_id" in primary_keys:
        assert table == "mcause_project_datasets", "project_id is not available for {table}"

    overlap = get_overlapping_nid_extracts(df, table, conn_def=conn_def, primary_keys=primary_keys)
    if len(overlap) > 0:
        if replace:
            if "project_id" in primary_keys:
                delete_project_nid_extracts_from_table(table, overlap, conn_def=conn_def)
            else:
                delete_nid_extracts_from_table(table, overlap, conn_def=conn_def)
        else:
            raise AssertionError(
                "These nid-extracts are already in {}: " "\n{}".format(table, overlap)
            )

    engine = ezfuncs.get_engine(conn_def=conn_def)
    conn = engine.connect()
    df.to_sql(table, conn, if_exists="append", index=False)
    conn.close()
    print("Successfully uploaded {} rows to {}".format(len(df), table))


def delete_phase_output(phase, nid, extract_type_id, sub_dirs=None, launch_set_id=None):
    if launch_set_id:
        assert launch_set_id == 0, "Should never delete anything in archive except launch set 0"

    filepath = construct_phase_path("FILEPATH")
    filepath = filepath[:-4]
    hdf = filepath + ".h5"
    csv = filepath + ".csv"
    parquet = filepath + ".parquet"
    filepaths = [csv, hdf, parquet]
    for filepath in filepaths:
        if os.path.exists(filepath):
            os.unlink(filepath)


def get_datestamp_modtime(filepath):
    return datetime.fromtimestamp(os.path.getmtime(filepath))


def get_launch_set_id_of_active(phase, nid, extract_type_id, sub_dirs=None):
    phase_path = construct_phase_path(
        phase, nid, extract_type_id, sub_dirs=sub_dirs, file_format="csv"
    )
    if not os.path.exists(phase_path):
        print(f"{phase_path} does not exist")
        return np.NaN
    curr_datestamp = get_datestamp_modtime(phase_path)

    data_folder = os.path.dirname(phase_path)
    arch_folder = os.path.join("FILEPATH")
    arch_files = os.listdir("FILEPATH")
    phase_arch_files = [f for f in arch_files if phase in f]

    close_time_launch_sets = []

    time_elapsed = timedelta(0, 3600)

    for phase_arch_file in phase_arch_files:
        test_datestamp = get_datestamp_modtime("FILEPATH")

        if abs(curr_datestamp - test_datestamp) < time_elapsed:
            pattern = "{phase}_([0-9_]+)\.csv".format(phase=phase)
            lsid_search = re.search(pattern, phase_arch_file)
            assert lsid_search is not None, "{} did not match {}".format(phase_arch_file, pattern)
            launch_set_id = lsid_search.group(1)
            close_time_launch_sets.append(launch_set_id)

    assert len(close_time_launch_sets) > 0, (
        "Cannot determine launch set id,"
        " the archive and current file"
        " timestamps are too far apart"
    )
    launch_set_id = max(close_time_launch_sets)
    return launch_set_id


def get_mcause_data(
    phase,
    nid=None,
    extract_type_id=None,
    project_id=None,
    source=None,
    nid_extract_records=None,
    location_id=None,
    sample=False,
    sample_num=None,
    force_rerun=False,
    data_type_id=None,
    year_id=None,
    code_system_id=None,
    iso3=None,
    region_id=None,
    launch_set_id=None,
    location_set_id=None,
    location_set_version_id=None,
    is_active=None,
    exec_function=None,
    exec_function_args=[],
    block_rerun=True,
    where_filter=None,
    sub_dirs=None,
    assert_all_available=False,
    verbose=False,
    attach_launch_set_id=False,
    refresh_id=None,
    usecols=None,
    cache_results=False,
):
    if attach_launch_set_id:
        assert launch_set_id is None, (
            "You want to pull data for a previous launch set and attach "
            "the launch set id for the current data; this information "
            "is incompatible."
        )

    if sample:
        assert sample_num is not None, (
            "If sampling datasets, must pass number of datasets to "
            "sample with sample_num argument. Got: {}".format(sample_num)
        )

    if launch_set_id:
        assert refresh_id is None, "Cannot specify both launch_set_id and refresh_id"

    dataset_filters = {
        "nid": nid,
        "extract_type_id": extract_type_id,
        "project_id": project_id,
        "nid_extract_records": nid_extract_records,
        "source": source,
        "location_id": location_id,
        "year_id": year_id,
        "data_type_id": data_type_id,
        "code_system_id": code_system_id,
        "iso3": iso3,
        "region_id": region_id,
        "location_set_id": location_set_id,
        "location_set_version_id": location_set_version_id,
        "is_active": is_active,
    }
    if verbose:
        print_log_message("Getting datasets to read")

    datasets = get_datasets(
        force_rerun=force_rerun,
        block_rerun=block_rerun,
        cache_results=cache_results,
        verbose=verbose,
        **dataset_filters,
    )
    if verbose:
        print_log_message("Got {} datasets".format(len(datasets)))

    pairs = datasets[["nid", "extract_type_id"]].drop_duplicates()

    if sample:
        if sample_num > len(pairs):
            warnings.warn(
                "Sample num of ({}) exceeded number of datasets "
                "({})".format(sample_num, len(pairs))
            )
            sample_num = len(pairs)
        pairs = pairs.sample(n=sample_num)
        if verbose:
            print_log_message("Restricted to {} random datasets".format(sample_num))
    nid_extract_pairs = [tuple(pair) for pair in list(pairs.to_records(index=False))]

    if verbose:
        print_log_message("Checking which datasets have available files")
    avail_nid_extracts = []
    for nid, extract_type_id in nid_extract_pairs:
        if check_output_exists(
            phase,
            nid,
            extract_type_id,
            sub_dirs=sub_dirs,
            refresh_id=refresh_id,
            project_id=project_id,
        ):
            avail_nid_extracts.append((nid, extract_type_id))

    bad_nid_extracts = list(set(nid_extract_pairs) - set(avail_nid_extracts))
    info_string = "Found {n} files to read data for.".format(n=len(avail_nid_extracts))
    if len(bad_nid_extracts) > 0:
        info_string = info_string + f" These nids were not available: \n{bad_nid_extracts}"
    if launch_set_id:
        info_string += " Using launch set id {}".format(launch_set_id)
    if assert_all_available and len(bad_nid_extracts) != 0:
        raise AssertionError(info_string)
    elif verbose:
        print_log_message(info_string)

    if len(avail_nid_extracts) == 0:
        raise AssertionError(
            "No files were available with given dataset " "filters: \n{}".format(dataset_filters)
        )
    if verbose:
        print_log_message(
            "Reading and appending {} data for {} "
            "nid-extracts".format(phase, len(avail_nid_extracts))
        )
    dfs = []
    if (usecols is not None) and verbose:
        print_log_message("Reading only a subset of the columns")
    for nid, extract_type_id in avail_nid_extracts:
        try:
            df = get_phase_output(
                phase,
                nid,
                extract_type_id,
                sub_dirs=sub_dirs,
                exec_function=exec_function,
                exec_function_args=exec_function_args,
                where_filter=where_filter,
                refresh_id=refresh_id,
                usecols=usecols,
                launch_set_id=launch_set_id,
                project_id=project_id,
            )

            check_cols = [x for x in df.columns if "multiple_cause_" in x]
            for col in check_cols:
                if set(df[col]) == set(["0000"]):
                    df = df.drop(col, axis=1)

            extra_filter_keys = ["location_id", "year_id"]
            for var in extra_filter_keys:
                if var in list(dataset_filters.keys()) and dataset_filters[var] is not None:
                    vals = dataset_filters[var]
                    if not isinstance(vals, list):
                        vals = [vals]
                    df = df.loc[df[var].isin(vals)]
                    if verbose:
                        print_log_message(
                            "Pruned dataset by {} filter, leaving "
                            "{} remaining".format(var, len(df))
                        )
        except KeyError as ke:
            fp_signature = f" [phase: {phase}, nid: {nid}, extract_type_id: {extract_type_id}]"
            etext = str(ke) + fp_signature
            if assert_all_available:
                raise KeyError(etext)
            else:
                warnings.warn(etext)
                continue
        if attach_launch_set_id:
            df["launch_set_id"] = get_launch_set_id_of_active(
                phase, nid, extract_type_id, sub_dirs=sub_dirs
            )
        dfs.append(df)
    print("Concatenating dataframes!")
    df = pd.concat(dfs, ignore_index=True, sort=True)

    if verbose:
        print_log_message("Constructed a dataset of {} rows".format(len(df)))

    return df


class McauseResult:

    processes = [
        "run_model",
        "predict",
        "remove_predictions_template",
        "calculate_counts",
        "age_loc_aggregation",
        "cause_aggregation",
        "calculate_props",
        "compile",
    ]
    process_suffix_dict = {
        "predictions_template": "_template",
        "predict": "",
        "calculate_counts": "_counts",
        "age_loc_aggregation": "_aggregate",
        "cause_aggregation": "_aggregate",
        "calculate_props": "_prop",
        "compile": "_summary",
        "squeeze": "_squeezed",
    }
    valid_intermediate_causes = list(MCoDMapper.int_cause_name_dict.keys()) + ["explicit_sepsis"]
    infectious_syndromes = MCoDMapper.infectious_syndromes
    assert not set(infectious_syndromes).issubset(
        valid_intermediate_causes
    ), "oops, syndromes and intermediate causes overlap, need another solution!"
    product_list = ["mortality", "rdp", "incidence", "attributable_burden"]

    def __init__(
        self,
        int_cause: str,
        end_product: str,
        process: str,
        year_id: int = None,
        cause_id=None,
        age_group_id: int = None,
        status: str = None,
        description: str = None,
        suffix: str = "parquet",
        conf: Configurator = None,
    ):
        self.process = process
        self.end_product = end_product
        self.int_cause = int_cause
        self.year_id = year_id
        self.cause_id = cause_id
        self.age_group_id = age_group_id
        self.status = status
        self.description = description
        self.suffix = suffix
        self.conf = conf or Configurator()
        assert self.conf.config_type in ["amr", "mcod"]
        if self.int_cause in self.infectious_syndromes:
            assert self.conf.config_type == "amr"
        assert self.process in self.processes + list(self.process_suffix_dict.keys())
        assert self.end_product in self.product_list
        self.construct_results_path()

    def get_model_version_path(self) -> PosixPath:
        assert (self.status in ["best", "latest"]) ^ (self.description is not None)
        if self.conf.config_type == "amr":
            if self.int_cause == "sepsis":
                model_step = "01_sepsis"
            else:
                model_step = "02_infectious_syndrome"
        else:
            model_step = ""
        base_dir = "FILEPATH"

        if self.status == "best":
            gbd_df = get_ids("gbd_round").assign(
                version_string=lambda x: "GBD" + x["gbd_round"].astype(str)
            )
            gbd_round = self.conf.get_id("gbd_round")
            if self.end_product == "rdp":
                gbd_round += 1
            best = gbd_df.loc[gbd_df["gbd_round_id"] == gbd_round, "version_string"]
            assert len(best) == 1, f"Found > 1 best version: {best}"
            assert (base_dir / best.iloc[0]).exists()
            return base_dir / best.iloc[0]
        elif self.status == "latest":
            version_dict = {}
            for model_dir in base_dir.iterdir():
                c_time = model_dir.stat().st_ctime
                version_dict.update({model_dir: c_time})
            most_recent = []
            for model_dir, c_time in version_dict.items():
                if (c_time == max(version_dict.values())) & ("GBD" not in model_dir.name):
                    most_recent.append(model_dir)
            assert len(most_recent) == 1, f"Found > 1 version: {most_recent}"
            return most_recent[0]
        else:
            return "FILEPATH"

    def construct_results_path(self):
        model_dir = "FILEPATH"
        self.parent_model_dir = "FILEPATH"
        if self.process in ["run_model", "predictions_template"]:
            if "by_age" in str(model_dir):
                model_dir = "FILEPATH"
        if self.process == "run_model":
            self.results_path = "FILEPATH"
        else:
            if (self.process == "compile") and (self.end_product == "rdp"):
                self.results_path = "FILEPATH"
            else:
                self.year_id = int(self.year_id)
                if self.cause_id != "_all_int_cause":
                    self.cause_id = int(self.cause_id)
                self.results_path = "FILEPATH"
            self.results_path = "FILEPATH"

    def make_results_directory(self):
        makedirs_safely("FILEPATH")

    def clear_results(self):
        if self.results_path.exists():
            self.results_path.unlink()

    def write_results(self, df: pd.DataFrame):
        self.make_results_directory()
        if self.results_path.suffix == ".h5":
            df.to_hdf("FILEPATH", key="model", format="fixed", mode="w")
        elif self.results_path.suffix == ".csv":
            df.to_csv("FILEPATH", index=False)
        elif self.results_path.suffix == ".parquet":
            df.to_parquet("FILEPATH", index=False)
        else:
            raise NotImplementedError

    def delete_last_step_results(self):

        current_path = "FILEPATH"

        if self.process == "age_loc_aggregation":
            file = current_path.replace("_aggregate", "_counts")
        elif self.process == "compile":
            file = current_path.replace("_summary", "_prop")
        else:
            file = []

        if file:
            if os.path.exists(file):
                os.unlink(file)

    def write_to_archive(self, df: pd.DataFrame):
        old_path = "FILEPATH"
        self.results_path = Path("FILEPATH")
        self.write_results(df)
        self.results_path = "FILEPATH"

    def read_results(self, usecols: List[str] = None) -> pd.DataFrame:
        assert self.results_path.suffix in [".h5", ".csv", ".parquet"]
        if self.results_path.suffix == ".h5":
            try:
                return pd.read_hdf("FILEPATH", key="model", columns=usecols)
            except OSError:
                print_log_message("H5 file not found, trying for historic csv...")
                return pd.read_csv("FILEPATH", usecols=usecols)
        elif self.results_path.suffix == ".csv":
            return pd.read_csv("FILEPATH", usecols=usecols)
        elif self.results_path.suffix == ".parquet":
            return pd.read_parquet("FILEPATH", columns=usecols)


def get_mcause_results(
    int_cause,
    end_product,
    year_id,
    cause_id,
    process,
    status=None,
    description=None,
    age_group_id=None,
    sex_id=None,
    location_id=None,
    exec_function=None,
    exec_function_args=None,
    usecols=None,
    suffix="parquet",
):
    exec_function_args = exec_function_args or {}
    year_id = wrap(year_id)
    cause_id = wrap(cause_id)

    dfs = []
    for year in year_id:
        for cause in cause_id:
            df = McauseResult(
                int_cause=int_cause,
                end_product=end_product,
                process=process,
                year_id=year,
                cause_id=cause,
                status=status,
                description=description,
                suffix=suffix,
                conf=CONF,
            ).read_results(usecols=usecols)
            if exec_function is not None:
                df = exec_function(df, **exec_function_args)
            dfs.append(df)
    df = pd.concat(dfs, sort=False)

    filters = {
        "age_group_id": age_group_id,
        "sex_id": sex_id,
        "location_id": location_id,
    }
    for var, values in filters.items():
        if values is not None:
            df = df.loc[df[var].isin(wrap(values))]
    return df


def setup_logger(log_title, log_path=None):
    logger = logging.getLogger("my_logger")

    logger.setLevel(logging.DEBUG)

    if log_path == None:
        user = getpass.getuser()
        log_path = "FILEPATH"

    fh = logging.FileHandler("FILEPATH")
    fh.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    fh.setFormatter(formatter)

    logger.addHandler(fh)

    return logger
