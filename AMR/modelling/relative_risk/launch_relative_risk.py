
import argparse
import os

import numpy as np
import pandas as pd

from cod_prep.claude.configurator import Configurator
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from mcod_prep.utils.mcause_io import makedirs_safely


class RelativeRiskCalculator(object):
    conf = Configurator("standard")

    def __init__(self, run_filters):
        assert (
                self.conf.config_type == "amr"
            ), "please change 'FILEPATH' to 'FILEPATH"

        self.processes = run_filters["processes"]
        self.description = run_filters["description"]
        self.sources = run_filters["sources"]
        self.custom_input_dir = run_filters["custom_input_dir"]
        self.isv = run_filters["isv"]
        self.oosv = run_filters["oosv"]
        self.burden = run_filters["burden"]
        self.log_base_dir = "FILEPATH"
        self.cor_num = 35


    def launch_data_prep(self): 
        worker = os.path.join("FILEPATH")
        jobname = f"rr_data_prep_{self.burden}_{self.description}"
        params = [self.description, self.burden] + self.sources 


        submit_mcod(
                jobname,
                "python",
                worker,
                cores=self.cor_num,
                verbose=True,
                logging=True,
                memory="150G",
                params=params,
                runtime="02:00:00",
            )


    def launch_mca(self, holds=[]): 
        worker = os.path.join("FILEPATH")
        jobname = f"rr_mca_{self.description}"
        params = [self.description]

        makedirs_safely("FILEPATH")

        submit_mcod(
                jobname,
                "r",
                worker,
                cores=1,
                verbose=True,
                logging=True,
                memory="100G",
                params=params,
                runtime="02:00:00",
                image="FILEPATH"
            )


    def launch_modeling(self, holds=[]): 
        if self.burden == "fatal":
            worker = os.path.join("FILEPATH")
            model_types = ["none", "all"]
        else: 
            worker = os.path.join("FILEPATH")
            model_types = ["none"]

        if self.custom_input_dir != "0": 
            input_dir = self.custom_input_dir
        else: 
            input_dir = self.description

        for syn_type in model_types: 

            if syn_type == "none": 
                save_dir = f"{self.description}_no_syndrome"
            else: 
                save_dir = self.description

            makedirs_safely("FILEPATH")
            makedirs_safely("FILEPATH")
            makedirs_safely("FILEPATH")
            makedirs_safely("FILEPATH")

            jobname = f"rr_modeling_{self.description}_{self.burden}_{syn_type}_syndrome"
            params = [input_dir, save_dir, self.isv, self.oosv, syn_type]

            submit_mcod(
                    jobname,
                    "r",
                    worker,
                    cores=1,
                    verbose=True,
                    logging=True,
                    memory="100G",
                    params=params,
                    runtime="02:00:00",
                )


    def launch(self):
        data_prep_jobs = []
        mca_jobs = []
        modeling_jobs = []

        if "data_prep" in self.processes:
            jid = self.launch_data_prep()
            data_prep_jobs.append(jid)
        if "mca" in self.processes:
            jid = self.launch_mca(holds=data_prep_jobs)
            mca_jobs.append(jid)
        if "modeling" in self.processes:
            jid = self.launch_modeling(holds=data_prep_jobs)
            modeling_jobs.append(jid)

    def check_output_exists(self):
        missing_processes = {}
        for process in self.processes:
            missing_sources = []
            if process == "data_prep":
                if self.burden == "fatal": 
                    sources = ["SOURCES"]
                    subfolder = "coefficients"
                    suffix = "coeff"
                else: 
                    sources = ["SOURCES"]
                    subfolder = "los"
                    suffix = "los"
                for source in sources: 
                    results_path = "FILEPATH"
                    if os.path.exists(results_path):
                        size = os.path.getsize(results_path)
                        if size <= 0:
                            missing_sources.append(source)
                            missing_processes.update({process: missing_sources})
                    else:
                        missing_sources.append(source)
                        missing_processes.update({process: missing_sources})
            elif process == "modeling":
                if self.burden == "fatal": 
                    model_types = [self.description, f"{self.description}_no_syndrome"]
                else: 
                    model_types = [f"{self.description}_no_syndrome"]
                for model_name in model_types: 
                    antibiotic_group = ['third_gen_ceph','carbapenem','fluoroquinolone','penicillin','aminopenicillin','beta_lactamase_inhibitor',
                    'anti_pseudomonal_penicillin','methicillin','vancomycin','fourth_gen_ceph','sulfa','aminoglycoside','macrolide']
                    for abx in antibiotic_group: 
                        results_path = "FILEPATH"
                        if os.path.exists(results_path):
                            size = os.path.getsize(results_path)
                            if size <= 0:
                                missing_sources.append(abx)
                                missing_processes.update({process: missing_sources})
                        else:
                            missing_sources.append(abx)
                            missing_processes.update({process: missing_sources})
        if len(missing_processes) > 0:
            print("Oops - The following were not available: {}".format(missing_processes))
        else:
            print("All sources/models are available!")

        return missing_processes


if __name__ == "__main__":
    processes = ["data_prep", "modeling", "mca"]
    parser = argparse.ArgumentParser(description="Launch relative risk analyses")
    
    parser.add_argument(
        "processes",
        help="processes to be run, e.g. data_prep, mca, modeling",
        type=str,
        nargs="+",
        choices=processes,
    )

    parser.add_argument(
        "description",
        help="a description of the model and how the data is split - ideally something like 'no_split_gram2' or 'age_location_split'",
        type=str,
    )

    parser.add_argument(
        "burden",
        help="Whether to run for fatal or nonfatal"
        'Default is to run for fatal.',
        type=str,
        default="fatal",
    )

    parser.add_argument(
        "--sources",
        help="Whether to run for specific sources"
        'Default is to run for all.',
        type=str,
        nargs="+",
        default=["all"],
    )

    parser.add_argument(
        "--custom_input_dir",
        help="If you want to use a specific version of data prep, you can specify the folder name of the data prep version you want",
        type=str,
        default="0"
    )

    parser.add_argument(
        "--isv",
        help="for modeling only: run in sample validation. Default is to not run."
        "Note: ISV not available for nonfatal yet",
        action="store_const",
        const=1,
        default=0,
    )

    parser.add_argument(
        "--oosv",
        help="for modeling only: run out of sample validation. Default is to not run."
        "Note: OOSV not available for nonfatal yet",
        action="store_const",
        const=1,
        default=0,
    )


    parser.add_argument(
        "--check", action="store_true", help="Do not launch, just check for outputs - available for data prep & modeling"
    )

    args = vars(parser.parse_args())
    check = args["check"]
    args.pop("check")
    launcher = RelativeRiskCalculator(args)
    print(f"You've submitted the following arguments: {args}")
    if check:
        launcher.check_output_exists()
    else:
        launcher.launch()
