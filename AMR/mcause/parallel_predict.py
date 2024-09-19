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
from mcod_prep.utils.mcause_io import setup_logger

CONF = Configurator()
NUM_WORKER = 38


def launch_predictionsworker(description, int_cause, parent_model_dir, year_id, cause_id):
    logger = setup_logger(f"{int_cause}_predict_{year_id}_{cause_id}")
    worker = "FILEPATH"

    logger.debug(f"Working on year {year_id} and cause {cause_id}")

    cmd = [
        "FILEPATH",
        "-s",
        worker,
        description,
        int_cause,
        parent_model_dir,
        year_id,
        cause_id,
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
                logger.info(f"R script output for year_id {year_id} and cause_id {cause_id}:{line}")
            for line in proc.stderr:
                logger.error(f"R script error for year_id {year_id} and cause_id {cause_id}:{line}")

            proc.wait()
            if proc.returncode != 0:
                logger.error(
                    f"R script failed with return code {proc.returncode} for year_id {year_id} and cause_id {cause_id}"
                )

    except Exception as e:
        logger.error(
            f"Exception occurred during subprocess call for year_id {year_id} and cause_id {cause_id}:\n{e}"
        )


def main(description, int_cause, parent_model_dir):
    array_df = pd.read_parquet("FILEPATH")
    print_log_message("Launching predictions by year and cause")

    with multiprocessing.Manager() as manager:
        with Pool(processes=NUM_WORKER) as pool:
            args_list = [
                (description, int_cause, parent_model_dir, str(row.year_id), str(row.cause_id))
                for row in array_df.itertuples(index=False)
            ]
            for _ in tqdm(pool.starmap(launch_predictionsworker, args_list)):
                pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format data and launch model.")
    parser.add_argument("description", type=str)
    parser.add_argument("int_cause", type=str)
    parser.add_argument("parent_model_dir", type=str)
    args = parser.parse_args()
    main(**vars(args))
