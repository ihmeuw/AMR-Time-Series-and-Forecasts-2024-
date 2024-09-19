import argparse
import itertools
import logging
import multiprocessing
import os
import shlex
import subprocess
import sys
from multiprocessing import Pool, current_process

import pandas as pd
from tqdm import tqdm

from cod_prep.claude.configurator import Configurator
from cod_prep.utils import print_log_message
from mcod_prep.calculate_counts import main_calculate_counts
from mcod_prep.utils.mcause_io import setup_logger

CONF = Configurator()
NUM_WORKER = 39


def launch_calculationsworker(description, int_cause, end_product, no_squeeze, year_id, cause_id):
    logger = setup_logger(f"calculate_counts_{int_cause}_{year_id}_{cause_id}")
    logger.debug(f"Working on year {year_id} and cause {cause_id}")

    # Load in the main of calculate counts
    main_calculate_counts(
        description, cause_id, year_id, end_product, int_cause, no_squeeze, logger
    )

    logger.debug(f"Finishing year {year_id} and cause {cause_id}")


def parallel_main(description, int_cause, end_product, year_id, no_squeeze, parent_model_dir):
    array_df = pd.read_parquet("FILEPATH")
    array_df = array_df.loc[array_df["year_id"] == year_id]
    print_log_message("Launching calculate counts by year and cause")

    with multiprocessing.Manager() as manager:
        with Pool(processes=NUM_WORKER) as pool:
            args_list = [
                (
                    description,
                    int_cause,
                    end_product,
                    no_squeeze,
                    row.year_id,
                    row.cause_id,
                )
                for row in array_df.itertuples(index=False)
            ]
            for _ in tqdm(pool.starmap(launch_calculationsworker, args_list)):
                pass


if __name__ == "__main__":
    description = str(sys.argv[1])
    int_cause = str(sys.argv[2])
    end_product = str(sys.argv[3])
    year_id = int(sys.argv[4])
    no_squeeze = str(sys.argv[5])
    parent_model_dir = str(sys.argv[6])

    assert no_squeeze in ["True", "False"]
    no_squeeze = no_squeeze == "True"

    print_log_message(f"No Squeeze: {no_squeeze}")
    parallel_main(description, int_cause, end_product, year_id, no_squeeze, parent_model_dir)
