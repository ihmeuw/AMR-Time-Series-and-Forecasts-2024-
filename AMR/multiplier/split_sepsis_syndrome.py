import pandas as pd
import logging
import sys
import os
import numpy as np
from cod_prep.claude.configurator import Configurator
from mcod_prep.calculate_counts import get_codcorrect_deaths
from mcod_prep.utils.mcause_io import get_mcause_results, setup_logger
from cod_prep.utils import print_log_message, report_if_merge_fail
from cod_prep.downloaders import (
    get_all_related_causes,
    get_current_cause_hierarchy
)
from amr_prep.utils.amr_io import AmrResult
from amr_prep.utils.misc import get_prepped_csbg_universe
CONF = Configurator()


class SepsisSyndromeSplitter(object):

    def __init__(self, year_id, cause_id, num_workers):
        self.dem_cols = ['location_id', 'sex_id', 'year_id', 'age_group_id']
        self.cache_kwargs = {'force_rerun': False, 'block_rerun': True, 'cache_results': False}
        self.conf = Configurator()
        self.draw_cols = ["draw_" + str(x) for x in range(0, self.conf.get_id('draw_num'))]
        self.cc_draw_cols = ["draw_" + str(x) + '_cc' for x in range(0, self.conf.get_id('draw_num'))]
        self.scl = pd.read_csv("FILEPATH")
        self.year_id = year_id
        self.cause_id = cause_id
        self.num_workers = num_workers

    def get_sepsis_by_syndrome(self):
        print_log_message("Reading sepsis draws")
        sepsis = get_mcause_results(
            int_cause='sepsis',
            end_product='mortality',
            year_id=self.year_id,
            cause_id=self.cause_id,
            process='squeeze',
            description=self.conf.get_id("sepsis_model"),
            usecols=self.dem_cols + self.draw_cols + ['point_estimate']
        )
        if self.scl.loc[self.scl['cause_id'] == self.cause_id, '100_sepsis'].item() == 1:
            sepsis[self.draw_cols + ['point_estimate']] = 1
         
        if self.cause_id in [322, 595]:
            print_log_message(
                'For GBD LRI and urinary_nephritis, no explicit syndrome outputs exist ' \
                'because they are simply the community syndrome definition'
            )
            if self.cause_id == 322:
                self.df = sepsis.assign(infectious_syndrome='L2_lower_respiratory_infection_comm')
            elif self.cause_id == 595:
                self.df = sepsis.assign(infectious_syndrome='L2_urinary_tract_infection_comm')
        elif self.cause_id in [394, 395, 396, 399]:
            self.df = sepsis.assign(infectious_syndrome='L2_sexually_transmitted_infection')
        else:
            print_log_message(f"Reading syndrome draws")
            best_models = pd.read_csv("FILEPATH")
            scll = pd.melt(
                self.scl,
                id_vars='cause_id',
                value_vars=best_models['int_cause'].unique().tolist(),
                var_name='int_cause',
                value_name='modeled'
            )
            int_causes = scll.loc[
                (scll['cause_id'] == cause_id) & (scll['modeled'] == 1),
                'int_cause'
            ].unique().tolist()
            best_models = best_models.loc[best_models['int_cause'].isin(int_causes), ]
            best_models = best_models.set_index('int_cause')['model_description'].to_dict()
            infsyn = pd.concat([
                get_mcause_results(
                    int_cause=syndrome,
                    end_product='mortality',
                    year_id=self.year_id,
                    cause_id=self.cause_id,
                    process='squeeze',
                    description=version,
                    usecols=self.dem_cols + self.draw_cols + ['point_estimate']
                ).assign(infectious_syndrome=syndrome)
                for syndrome, version in best_models.items()
            ], sort=False)

            print_log_message("Multiplying sepsis by syndrome")
            assert np.allclose(infsyn.groupby(['year_id', 'location_id', 'age_group_id', 'sex_id'])
                               [self.draw_cols + ['point_estimate']].sum(), 1)
            self.df = pd.merge(
                sepsis, infsyn, on=self.dem_cols,
                validate='one_to_many', how='outer', indicator=True,
                suffixes=('_sep', '_syn')
            )
            assert (self.df._merge == 'both').all()
            self.df[self.draw_cols] = pd.DataFrame(
                self.df.filter(regex=r'draw_\d{1,3}_sep').to_numpy()
                * self.df.filter(regex=r'draw_\d{1,3}_syn').to_numpy(),
                index=self.df.index
            )
            self.df['point_estimate'] = self.df['point_estimate_sep'] * self.df['point_estimate_syn']
        
        self.df = self.df[self.dem_cols + ['infectious_syndrome'] + self.draw_cols + ['point_estimate']]
        assert self.df.notnull().values.all()
        assert (self.df[self.draw_cols] <= 1).values.all()

    def multiply_codcorrect(self):
        print_log_message(f"Reading CodCorrect draws")
        ages = self.df.age_group_id.unique().tolist()
        locations = self.df.location_id.unique().tolist()
        sexes = self.df.sex_id.unique().tolist()
        logger = setup_logger('split_sepsis_syndrome_' + str(self.year_id) + '_' + str(self.cause_id))
        print(self.num_workers)
        cc = get_codcorrect_deaths(
            self.cause_id, self.year_id, ages,
            locations, sexes, logger, num_workers=self.num_workers
        )
        cc = cc[self.dem_cols + self.draw_cols + ['point_estimate']]

        print_log_message("Multiplying in CodCorrect deaths")
        self.df = self.df.merge(
            cc, how='left', on=self.dem_cols, indicator=True,
            suffixes=('', '_cc'), validate='many_to_one'
        )
        if (self.conf.get_id("codcorrect_version") == 135) & (self.cause_id in [656, 657, 668]):
            self.df = self.df.loc[~self.df.age_group_id.isin([2, 3])]
        assert (self.df._merge == 'both').all()
        self.df[self.draw_cols] = pd.DataFrame(
            self.df[self.draw_cols].to_numpy()
            * self.df[self.cc_draw_cols].to_numpy(),
            index=self.df.index
        )
        self.df['point_estimate'] = self.df['point_estimate'] * self.df['point_estimate_cc']
        self.df = self.df[self.dem_cols + ['infectious_syndrome'] + self.draw_cols + ['point_estimate']]

    def run_split(self):
        print_log_message(f"Running split for cause {self.cause_id}, year {self.year_id}")
        self.get_sepsis_by_syndrome()
        self.multiply_codcorrect()

        if 'L2_mix_respiratory_infection' in self.df['infectious_syndrome'].unique().tolist():
            self.df.loc[
                self.df['infectious_syndrome'] == 'L2_mix_respiratory_infection',
                'infectious_syndrome'
            ] = 'L2_lower_respiratory_infection_hosp'
            self.df = self.df.groupby(
                self.dem_cols + ['infectious_syndrome'], as_index=False
            )[self.draw_cols + ['point_estimate']].sum()

        self.df['cause_id'] = self.cause_id
        self.df['hosp'] = 'all'
        self.df.loc[self.df.infectious_syndrome.str.endswith("_hosp"), 'hosp'] = 'hospital'
        self.df.loc[self.df.infectious_syndrome.str.endswith("_comm"), 'hosp'] = 'community'
        self.df['infectious_syndrome'] = self.df.infectious_syndrome.str.replace('_hosp', '')
        self.df['infectious_syndrome'] = self.df.infectious_syndrome.str.replace('_comm', '')

        assert self.df.notnull().values.all()

    def save_output(self):
        print_log_message("Saving outputs")
        for syndrome, syndrome_df in self.df.groupby('infectious_syndrome'):
            print_log_message(f"Working on syndrome {syndrome}")
            AmrResult(
                process='split_sepsis_syndrome',
                burden_type='fatal',
                year_id=self.year_id,
                cause_id=self.cause_id,
                infectious_syndrome=syndrome
            ).write_results(syndrome_df)


if __name__ == '__main__':
    num_workers = int(sys.argv[1])
    task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        base_dir = "FILEPATH"
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        cause_id = int(task_row['cause_id'])
    else:
        year_id = int(sys.argv[2])
        cause_id = int(sys.argv[3])

    splitter = SepsisSyndromeSplitter(
        year_id=year_id, cause_id=cause_id, num_workers=num_workers
    )
    splitter.run_split()
    splitter.save_output()
