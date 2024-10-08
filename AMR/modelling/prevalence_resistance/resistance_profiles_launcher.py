import pandas as pd
from cod_prep.utils import (
    report_duplicates, wait_for_job_ids,
    print_log_message
)
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from cod_prep.claude.configurator import Configurator
import getpass
from cod_prep.downloaders import get_current_location_hierarchy
user = getpass.getuser()

CONF = Configurator()

IN_DIR = "FILEPATH"

corrs = pd.read_csv(
    f'{IN_DIR}/subset_pearson_correlations.csv')
corrs_salmonella = pd.read_csv(
    f'{IN_DIR}/salmonella_pearson_corrs.csv')
corrs = corrs.append(corrs_salmonella)
corrs_staph = pd.read_csv(
    f'{IN_DIR}/staph_pearson_corrs.csv')
corrs = corrs.append(corrs_staph)
studycombos = pd.read_csv('FILEPATH'.format(user), encoding = 'latin-1')
bugs = [ele for ele in reversed(corrs['pathogen'].value_counts().index.tolist())]
bugs = ['staphylococcus_aureus']

locs = get_current_location_hierarchy(
    location_set_version_id=CONF.get_id('location_set_version'),
)

countries = list(locs.loc[locs['level'] == 3, 'location_id'])

log_base_dir = "FILEPATH"
for bug in bugs:
    print_log_message(f"Submitting jobs for {bug} models by location...")
    for loc in countries:
        worker = "FILEPATH"
        params = [bug, loc]
        jobname = f"model_{bug}_{loc}"
        submit_mcod(
            jobname, language='python', worker=worker,
            cores=1, memory="10G", params=params,
            runtime="10:00:00", logging=True,
            log_base_dir=log_base_dir
        )
    print_log_message("Finished submitting!")
