import argparse
import getpass
from builtins import object

import pandas as pd

from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import (
    get_ages,
    get_cause_map,
    get_current_cause_hierarchy,
    get_current_location_hierarchy,
    get_map_version,
    get_remove_decimal,
)
from cod_prep.utils import cod_timestamp, print_log_message
from mcod_prep.mcod_mapping import MCoDMapper
from mcod_prep.utils.mcause_io import check_output_exists, delete_phase_output, makedirs_safely
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from mcod_prep.utils.nids import add_nid_metadata, get_datasets


class MCauseLauncher(object):
    conf = Configurator()
    assert conf.config_type in ["mcod", "amr", "mcod_for_covid", "mcod_rdp"]
    cache_options = {
        "force_rerun": True,
        "block_rerun": False,
        "cache_dir": "standard",
        "cache_results": True,
        "verbose": True,
    }
    source_memory_dict = {
        "TWN_MOH": "12G",
        "MEX_INEGI": "100G",
        "BRA_SIM": "80G",
        "USA_NVSS": "80G",
        "ITA_ISTAT": "80G",
        "COL_DANE": "10G",
        "ITA_FVG_LINKAGE": "10G",
        "NZL_LINKAGE_NMDS": "15G",
        "ZAF_STATSSA": "50G",
        "NZL_NMDS": "350G",
        "MNG_H_INFO": "80G",
        "BRA_SIH": "80G",
        "AUT_HDD": "300G",
        "KGZ_MHIF": "3G",
        "ITA_HID": "150G",
        "USA_SID": "150G",
        "MEX_SAEH": "80G",
        "CAN_DAD": "1G",
        "USA_NHDS": "80G",
        "fl_doh": "15G",
    }

    def __init__(self, run_filters):
        self.run_filters = run_filters
        self.location_set_version_id = self.conf.get_id("location_set_version")
        self.cause_set_version_id = self.conf.get_id("cause_set_version")
        self.mcod_code = "FILEPATH"
        self.datasets = self.prep_run_filters()
        self.log_base_dir = "FILEPATH"

    def prep_run_filters(self):

        datasets_kwargs = {
            "force_rerun": True,
            "block_rerun": False,
            "cache_results": True,
            "is_active": True,
            "project_id": self.conf.get_id("project"),
        }
        datasets_kwargs.update(
            {
                k: v
                for k, v in self.run_filters.items()
                if k not in ["intermediate_causes", "phases"]
            }
        )
        index_cols = ["nid", "extract_type_id"]
        datasets = (
            get_datasets(**datasets_kwargs)
            .drop_duplicates(index_cols)
            .set_index(index_cols)[["year_id", "code_system_id", "source", "data_type_id"]]
        )
        gbdr_args = {}
        if self.conf.config_type == "amr":
            gbdr_args = {"gbd_round_id": self.conf.get_id("most_current_gbd_round")}
        datasets = datasets.assign(
            code_map_version_id=datasets["code_system_id"].apply(
                lambda x: get_map_version(x, "YLL", "best", **gbdr_args)
            ),
            remove_decimal=datasets["code_system_id"].apply(lambda x: get_remove_decimal(x)),
        )
        return datasets

    def cache_resources(self, cache_functions_to_run_with_args):
        for cache_function, kwargs in cache_functions_to_run_with_args:
            function_name = cache_function.__name__
            cache_exists = cache_function(only_check_cache_exists=True, verbose=True, **kwargs)
            if cache_exists:
                print(f"No need to recache method {function_name} with args: {kwargs}")
            else:
                print(f"Running {function_name} with args: {kwargs}")
                kwargs.update(self.cache_options)
                cache_function(**kwargs)

    def launch_format_map(
        self,
        year,
        source,
        int_cause,
        code_system_id,
        code_map_version_id,
        nid,
        extract_type_id,
        data_type_id,
    ):
        delete_phase_output("format_map", nid, extract_type_id, sub_dirs=int_cause)
        worker = "FILEPATH"
        params = [
            int(year),
            source,
            int_cause,
            int(code_system_id),
            int(code_map_version_id),
            int(nid),
            int(extract_type_id),
            int(data_type_id),
        ]
        jobname = f"mcause_format_map_{source}_{year}_{int_cause}"
        try:
            memory = self.source_memory_dict[source]
        except KeyError:
            print(f"{source} is not in source_memory_dict. Trying with 50G.")
            memory = "50G"
        if source in ["USA_SID", "AUT_HDD", "USA_NVSS", "ITA_HID", "BRA_SIH", "NZL_NMDS"]:
            runtime = "24:00:00"
        else:
            runtime = "02:00:00"
        if source in ["BRA_SIM", "COL_DANE", "ITA_ISTAT", "MEX_INEGI", "ZAF_STATSSA"]:
            jdrive = True
        else:
            jdrive = False
        jid = submit_mcod(
            jobname,
            "python",
            worker,
            cores=1,
            memory=memory,
            params=params,
            verbose=True,
            queue="long.q",
            logging=True,
            jdrive=jdrive,
            runtime=runtime,
            log_base_dir=self.log_base_dir,
        )
        return jid

    def launch_redistribution(
        self,
        nid,
        extract_type_id,
        code_system_id,
        code_map_version_id,
        remove_decimal,
        data_type_id,
        int_cause,
        holds=[],
        source=None,
    ):
        delete_phase_output("redistribution", nid, extract_type_id, sub_dirs=int_cause)
        worker = "FILEPATH"
        jobname = f"{int_cause}_redistribution_{nid}_{extract_type_id}"
        params = [
            nid,
            extract_type_id,
            self.cause_set_version_id,
            self.location_set_version_id,
            code_system_id,
            code_map_version_id,
            remove_decimal,
            int(data_type_id),
            int_cause,
        ]

        if int_cause in ["cross_tabulate", "standard_map"]:
            rd_source_mem_dict = {
                "ZAF_STATSSA": "150G",
                "NZL_NMDS": "300G",
                "BRA_SIH": "300G",
                "MNG_H_INFO": "300G",
                "AUT_HDD": "150G",
                "MEX_SAEH": "150G",
                "MEX_INEGI": "400G",
                "USA_SID": "100G",
                "USA_NHDS": "12G",
                "USA_NVSS": "400G",
                "BRA_SIM": "400G",
            }

            try:
                memory = rd_source_mem_dict[source]
            except KeyError:
                print(f"{source} is not in rd_source_mem_dict. Trying with 50G.")
                memory = "50G"
        else:
            rd_big_sources = [
                "USA_NVSS",
                "NZL_NMDS",
                "BRA_SIM",
                "BRA_SIH",
                "MNG_H_INFO",
                "AUT_HDD",
                "ITA_HID",
                "USA_SID",
                "MEX_SAEH",
                "USA_NHDS",
            ]
            if source is not None:
                if source in rd_big_sources:
                    memory = "300G"
                else:
                    memory = "50G"

        print(f"Running {source} with {memory}")
        runtime = "24:00:00"

        submit_mcod(
            jobname,
            "python",
            worker,
            cores=1,
            memory=memory,
            params=params,
            holds=holds,
            verbose=True,
            queue="long.q",
            logging=True,
            runtime=runtime,
            log_base_dir=self.log_base_dir,
        )

    def launch_reassign_injuries(self, int_cause, source):
        worker = "FILEPATH"
        if len(source) == 8:
            memory = "40G"
            cores = 17
        else:
            memory = "20G"
            cores = 10
        params = ["".join(int_cause), " ".join(source)]
        job_name = f"mcause_reassign_injuries_{'_'.join(source)}_{''.join(int_cause)}"
        submit_mcod(
            job_name,
            "python",
            worker,
            cores=cores,
            memory=memory,
            params=params,
            logging=True,
            jdrive=True,
            verbose=True,
            log_base_dir=self.log_base_dir,
        )

    def launch(self):
        conf = Configurator()
        if self.run_filters["intermediate_causes"] == ["all_syndromes"]:
            self.run_filters["intermediate_causes"] = conf.get_id("inf_syns")
        cache_functions_to_run_with_args = [
            (get_current_cause_hierarchy, {"cause_set_version_id": self.cause_set_version_id}),
            (get_ages, {}),
            (
                get_current_location_hierarchy,
                {"location_set_version_id": self.location_set_version_id},
            ),
        ]
        for code_map_version_id in list(self.datasets.code_map_version_id.unique()):
            cache_functions_to_run_with_args.append(
                (get_cause_map, {"code_map_version_id": code_map_version_id})
            )
        self.cache_resources(cache_functions_to_run_with_args)

        makedirs_safely(self.log_base_dir)

        if self.conf.get_id("project") == 2:
            allowed_int_causes = MCoDMapper.infectious_syndromes + ["sepsis"]
        elif self.conf.get_id("project") == 4:
            allowed_int_causes = list(MCoDMapper.int_cause_name_dict.keys()) + MCoDMapper.inj_causes
        elif self.conf.get_id("project") == 3:
            allowed_int_causes = ["lri_corona", "cross_tabulate"]
        elif self.conf.get_id("project") == 1:
            allowed_int_causes = ["cross_tabulate", "pathogen"]

        not_allowed = set(self.run_filters["intermediate_causes"]) - set(allowed_int_causes)
        assert set(self.run_filters["intermediate_causes"]).issubset(
            allowed_int_causes
        ), f"The following int_causes are not allowed in {self.conf.config_type}: \
            {not_allowed}"

        format_map_jobs = []
        for row in self.datasets.itertuples():
            nid, extract_type_id = row.Index
            for int_cause in self.run_filters["intermediate_causes"]:
                if "format_map" in self.run_filters["phases"]:
                    jid = self.launch_format_map(
                        row.year_id,
                        row.source,
                        int_cause,
                        row.code_system_id,
                        row.code_map_version_id,
                        nid,
                        extract_type_id,
                        row.data_type_id,
                    )
                    format_map_jobs.append(jid)
                if "redistribution" in self.run_filters["phases"]:
                    self.launch_redistribution(
                        nid,
                        extract_type_id,
                        row.code_system_id,
                        row.code_map_version_id,
                        row.remove_decimal,
                        row.data_type_id,
                        int_cause,
                        holds=format_map_jobs,
                        source=row.source,
                    )
        if "reassign_injuries" in self.run_filters["phases"]:
            for int_cause in self.run_filters["intermediate_causes"]:
                self.launch_reassign_injuries(int_cause, list(self.datasets.source.unique()))

    def check_output_exists(self):
        conf = Configurator()
        if self.run_filters["intermediate_causes"] == ["all_syndromes"]:
            self.run_filters["intermediate_causes"] = conf.get_id("inf_syns")
        failures = pd.DataFrame()
        for row in self.datasets.itertuples():
            nid, extract_type_id = row.Index
            for int_cause in self.run_filters["intermediate_causes"]:
                for phase in self.run_filters["phases"]:
                    if not check_output_exists(phase, nid, extract_type_id, sub_dirs=int_cause):
                        failures = failures.append(
                            {
                                "nid": nid,
                                "extract_type_id": extract_type_id,
                                "int_cause": int_cause,
                                "phase": phase,
                            },
                            ignore_index=True,
                        )
        if len(failures) == 0:
            print_log_message("Outputs all exist!")
        else:
            failures = add_nid_metadata(failures, ["source"])
            print_log_message(f"The following nid/etid/phase/int causes failed: \n {failures}")
            user = getpass.getuser()
            failures.to_csv("FILEPATH", index=False)

        return failures


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="NID/etids for which to run formatting, mapping, and redistribution"
    )
    parser.add_argument(
        "intermediate_causes",
        help="intermediate cause(s) of interest",
        nargs="+",
        choices=MCoDMapper.possible_int_causes + ["all_syndromes"],
        type=str,
    )
    parser.add_argument(
        "phases",
        help="data processing phases",
        type=str,
        nargs="+",
        choices=["format_map", "redistribution", "reassign_injuries"],
    )
    parser.add_argument("--iso3", nargs="*")
    parser.add_argument("--code_system_id", nargs="*", type=int)
    parser.add_argument("--data_type_id", nargs="*", type=int)
    parser.add_argument("--year_id", type=int, nargs="*")
    parser.add_argument("--source", nargs="*")
    parser.add_argument("--nid", nargs="*", type=int)
    parser.add_argument("--check", action="store_true")
    args = vars(parser.parse_args())
    check = args["check"]
    args.pop("check")
    launcher = MCauseLauncher(args)
    print(f"You've submitted the following arguments: {args}")
    if check:
        launcher.check_output_exists()
    else:
        launcher.launch()
