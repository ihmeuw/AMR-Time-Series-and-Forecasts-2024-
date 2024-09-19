import argparse
import itertools
import logging
import multiprocessing
import os
import shlex
import subprocess
from multiprocessing import Pool, current_process

import pandas as pd
from tqdm import tqdm

from cod_prep.claude.configurator import Configurator
from cod_prep.utils import print_log_message
from mcod_prep.utils.mcause_io import setup_logger, makedirs_safely

CONF = Configurator()
NUM_WORKER = 33

ucod_sources = ["SOURCES"] 
other_sources = ["SOURCES"]
biggest_source = ["SOURCES"]
fatal_sources = ucod_sources + other_sources + biggest_source 

nidmeta = pd.read_csv("SOURCES")

nonfatal_sources = nidmeta.loc[(nidmeta["active_5"] == 1) 
& (nidmeta["source"].notnull()), "source"].unique().tolist() + ["SOURCE"]

def launch_data_prep_worker(source, description, burden):

    logger = setup_logger(f"{source}_{description}_{burden}_data_prep")
    logger.debug(f"Working on {source}")

    if source != "SOURCE": 
        worker = "FILEPATH"
        cmd = [
            "FILEPATH",
            "-s",
            worker,
            source, 
            description,
            burden,
        ]
    else: 
        worker = "FILEPATH"
        cmd = [
            "FILEPATH",
            "-s",
            worker,
            description,
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
                logger.info(f"R script output for {source}:{line}")
            for line in proc.stderr:
                logger.error(f"R script error for {source}:{line}")

            proc.wait()
            if proc.returncode != 0:
                logger.error(
                    f"R script failed with return code {proc.returncode} for {source}"
                )

    except Exception as e:
        logger.error(
            f"Exception occurred during subprocess call for {source}:\n{e}"
        )


def main(description, burden, sources):

    if (sources == ["all"]) & (burden == "fatal"):  
        sources = fatal_sources
    elif (sources == ["all"]) & (burden == "nonfatal"): 
        sources = nonfatal_sources 

    makedirs_safely("FILEPATH")
    makedirs_safely("FILEPATH")
    
    print_log_message("Launching data prep by source")

    with multiprocessing.Manager() as manager:
        with Pool(processes=NUM_WORKER) as pool:
            args_list = [
                (source, description, burden)
                for source in sources
            ]
            for _ in tqdm(pool.starmap(launch_data_prep_worker, args_list)):
                pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format data and launch model.")
    parser.add_argument("description", type=str)
    parser.add_argument("burden", type=str)
    parser.add_argument("sources", nargs='+')
    args = parser.parse_args()
    main(**vars(args))