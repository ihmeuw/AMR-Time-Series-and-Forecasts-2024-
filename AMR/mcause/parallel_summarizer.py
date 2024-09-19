import multiprocessing
import sys
from multiprocessing import Pool, current_process

import pandas as pd
from tqdm import tqdm

from cod_prep.claude.configurator import Configurator
from cod_prep.utils import print_log_message
from mcod_prep.age_loc_aggregation import main_age_loc
from mcod_prep.back_calculate_proportions import main_props
from mcod_prep.cause_aggregation import main_cause_agg
from mcod_prep.compile_burden import main_compile
from mcod_prep.utils.mcause_io import makedirs_safely, setup_logger

CONF = Configurator()
NUM_WORKER = 39


def launch_summarizeworker(
    process, int_cause, description, end_product, cause_id, year_id, parent_model_dir
):
    logger = setup_logger(f"{int_cause}_{process}_{year_id}_{cause_id}")
    logger.debug(f"Working on {year_id}, {cause_id}")

    if process == "cause_aggregation":
        logger.debug("Aggregating Causes")
        main_cause_agg(description, year_id, int_cause, end_product, cause_id, logger)
    elif process == "age_loc_aggregation":
        logger.debug("Aggregating ages & locations")
        main_age_loc(int_cause, description, end_product, cause_id, year_id, logger)
    elif process == "back_calculate_props":
        logger.debug("Back-calculating proportions")
        makedirs_safely(f"{parent_model_dir}/bad_props")
        main_props(description, cause_id, year_id, end_product, int_cause, parent_model_dir, logger)
    else:
        logger.debug("Compiling")
        main_compile(description, end_product, int_cause, year_id, cause_id, logger)

    logger.debug("Finished!")


def main(process, description, end_product, int_cause, parent_model_dir, level, year):
    print_log_message(f"Launching {process} by year and cause")
    if process == "cause_aggregation":
        array_df = pd.read_parquet("FILEPATH")
    elif process == "age_loc_aggregation":
        array_df = pd.read_parquet("FILEPATH")
    elif process == "back_calculate_props":
        array_df = pd.read_parquet("FILEPATH")
        array_df = array_df.loc[array_df["year_id"] == year]
    else:
        array_df = pd.read_parquet("FILEPATH")

    with multiprocessing.Manager() as manager:
        with Pool(processes=NUM_WORKER) as pool:
            args_list = [
                (
                    process,
                    int_cause,
                    description,
                    end_product,
                    row.cause_id,
                    row.year_id,
                    parent_model_dir,
                )
                for row in array_df.itertuples(index=False)
            ]
            for _ in tqdm(pool.starmap(launch_summarizeworker, args_list)):
                pass


if __name__ == "__main__":
    process = str(sys.argv[1])
    description = str(sys.argv[2])
    end_product = str(sys.argv[3])
    int_cause = str(sys.argv[4])
    parent_model_dir = str(sys.argv[5])

    if process == "cause_aggregation":
        level = int(sys.argv[6])
    else:
        level = None

    if process == "back_calculate_props":
        year = int(sys.argv[6])
    else:
        year = None

    main(process, description, end_product, int_cause, parent_model_dir, level, year)
