import pandas as pd
import argparse
from cod_prep.utils import (
    report_duplicates, wait_for_job_ids,
    print_log_message
)
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from cod_prep.claude.configurator import Configurator
from modelling.pathogen.run_model import PathogenNetwork

CONF = Configurator()
HOLDS = []

def main(args):
    config_file = pd.read_excel(f"{FILEPATH}/config.xlsx", sheet_name='run')
    report_duplicates(config_file, ['model_version', 'infectious_syndrome'])
    worker = f"{FILEPATH}/run_model.py"
    print_log_message(f"Submitting jobs for {len(config_file)} models...")
    for model in config_file.to_dict('rows'):
        params = [
            model['model_version'],
            model['infectious_syndrome'],
            args.read_data_cache,
            args.read_model_cache,
            args.out_of_sample_validation
        ]
        jobname = f"model_{model['model_version']}_{model['infectious_syndrome']}"
        if model['infectious_syndrome'] == 'L2_lower_respiratory_infection':
            mem = "520G"
        else:
            mem = "250G"
        submit_mcod(
            jobname, language='python', worker=worker,
            cores=10, memory=mem, params=params,
            runtime="4:30:00", logging=True,
            log_base_dir=str(
                PathogenNetwork.out_dir / model['infectious_syndrome'] /
                model['model_version']
            ),
            holds=HOLDS
        )
    print_log_message("Finished submitting!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Launch pathogen network model')
    parser.add_argument('-rdc', '--read_data_cache', action='store_true')
    parser.add_argument('-rmc', '--read_model_cache', action='store_true')
    parser.add_argument('-oosv', '--out_of_sample_validation', action='store_true')
    args = parser.parse_args()
    main(args)
