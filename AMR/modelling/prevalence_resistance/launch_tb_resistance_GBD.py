import sys
import os
import pandas as pd
import glob
import numpy as np 
from cod_prep.downloaders.locations import add_location_metadata
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from amr_prep.utils.misc import get_prepped_csbg_universe
from cod_prep.claude.configurator import Configurator

CONF = Configurator('standard')

if __name__ == "__main__":
    for year in list(range(1990, 2022)):
        worker = CONF.get_directory('amr_repo') + 'filepath/prep_tb_resistance_GBD.py'
        jobname = f"prep_tb_GBD_prev_res_" + str(year)
        print('submitting ' + str(year))
        submit_mcod(
            jobname, language='python', worker=worker,
            cores=5, memory="5G",
            runtime='00:30:00',
            params=[year], logging=True, queue='all.q',
        )
