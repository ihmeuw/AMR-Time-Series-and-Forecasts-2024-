import pandas as pd
import sys
import os
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import AmrResult
from cod_prep.utils import (
    print_log_message,
    report_if_merge_fail,
    report_duplicates
)
from cod_prep.downloaders.ages import add_age_metadata
from cod_prep.downloaders.locations import add_location_metadata, get_current_location_hierarchy
from amr_prep.utils.misc import get_pred_ex, get_prepped_csbg_universe
import numpy as np
CONF = Configurator()


class AmrBurdenCalculator(object):

    def __init__(self, burden_type, year_id, cause_id, infectious_syndrome,
                 pathogen):
        self.burden_type = burden_type
        self.year_id = year_id
        self.cause_id = cause_id
        self.infectious_syndrome = infectious_syndrome
        self.pathogen = pathogen
        self.conf = Configurator()
        self.cache_kwargs = {'force_rerun': False, 'block_rerun': True, 'cache_results': False}
        self.dem_cols = ['location_id', 'age_group_id', 'sex_id', 'year_id']
        self.draw_cols = ["draw_" + str(x) for x in range(0, self.conf.get_id('draw_num'))]
        self.csbg = get_prepped_csbg_universe(burden_type=self.burden_type, **self.cache_kwargs)

    def calculate_ylls(self, df):
        df = df.copy()
        pred_ex = get_pred_ex(codcorrect_version=self.conf.get_id("codcorrect_version"))
        pred_ex_id_cols = ['location_id', 'age_group_id', 'sex_id', 'year_id']
        df = df.merge(
            pred_ex, how='left', on=pred_ex_id_cols,
            validate='many_to_one'
        )
        report_if_merge_fail(df, 'pred_ex', pred_ex_id_cols)
        draw_cols = self.draw_cols
        df[draw_cols] = pd.DataFrame(
            df.filter(regex='draw_\\d{1,3}$').to_numpy()
            * df[['pred_ex']].to_numpy(),
            index=df.index
        )
        df['point_estimate'] = df['point_estimate'] * df['pred_ex']
        df = df.drop('pred_ex', axis='columns')
        return df

    def calculate_ylds(self, df, ypc):
        df = df.copy()
        merge_cols = [
            'location_id', 'year_id', 'age_group_id', 'sex_id',
        ]
        if 'iso3' not in df.columns:
            df = add_location_metadata(df, 'iso3', **self.cache_kwargs)
            lh = get_current_location_hierarchy(**self.cache_kwargs)
            iso3_map = lh.loc[lh['level'] == 3, ['iso3', 'location_id']].drop_duplicates()
            iso3_map.rename(columns={'location_id': 'country_loc_id'}, inplace=True)
            df = df.merge(iso3_map, on=['iso3'], how='left', validate='many_to_one')
            assert df['country_loc_id'].notnull().values.all()
            ypc.rename(columns={'location_id': 'country_loc_id'}, inplace=True)
            merge_cols.remove('location_id')
            merge_cols.append('country_loc_id')
        df = df.merge(
            ypc, how='left', on=merge_cols, validate='many_to_one',
            indicator=True, suffixes=('', '_ypc')
        )
        assert (df._merge == 'both').all()
        df[self.draw_cols] = pd.DataFrame(
            df[self.draw_cols].to_numpy() * df.filter(regex='draw_\\d{1,3}_ypc').to_numpy(),
            index=df.index
        )
        df['point_estimate'] = df['point_estimate'] * df['pe_ylds_case']
        assert df[self.draw_cols].notnull().values.all()
        df = df.drop([f'draw_{x}_ypc' for x in range(0, CONF.get_id('draw_num'))] + ['_merge'], axis='columns')
        return df

    def multiply_prop(self, df, prop):
        assert (prop.age_group_id == 22).all()
        df['agg_age_group_id'] = 22
        df['parent_sex_id'] = 3
        prop = prop.rename(columns={'age_group_id': 'agg_age_group_id', 'sex_id': 'parent_sex_id'})
        prop_merge_cols = [
            'agg_age_group_id', 'parent_sex_id', 'location_id', 'year_id', 'pathogen',
            'measure_id'
        ]

        report_duplicates(prop, prop_merge_cols + ['counterfactual', 'abx_set', 'abx_class'])
        dp_prop = prop.loc[prop['abx_set'] == 'all', ]
        ndp_prop = prop.loc[prop['abx_set'] != 'all', ] 
        dfs = []
        if self.burden_type == 'fatal':
            measure = 1
        elif self.burden_type == 'nonfatal':
            measure = 6
        before_check = df.loc[
            (df['measure_id'] == measure),
            'point_estimate'
        ].sum()
        for prop in [dp_prop, ndp_prop]:
            if len(prop) == 0:
                continue
            prop = add_location_metadata(prop, 'level', **self.cache_kwargs)
            dfc = df.copy()
            if (prop['level'] == 3).values.all():
                dfc['detailed_location'] = dfc['location_id']
                dfc = add_location_metadata(dfc, ['path_to_top_parent'], **self.cache_kwargs)
                dfc['country_id'] = dfc['path_to_top_parent'].str.split(',', expand=True)[3]
                dfc['location_id'] = dfc['country_id'].astype(int)
                dfc.drop(columns=['path_to_top_parent', 'country_id'], inplace=True)
            prop.drop(columns=['level'], inplace=True)
            assert set(prop.location_id) == set(dfc.location_id)
            dfc = dfc.set_index(sorted([c for c in dfc if ('draw' not in c) & (c != 'point_estimate')]))
            prop = prop.set_index(sorted([c for c in prop if ('draw' not in c) & (c != 'point_estimate')]))
            dfc = dfc.mul(prop).reset_index(drop=False)
            dfs.append(dfc)
        df = pd.concat(dfs)
        df = df.reset_index(drop=True)

        after_check = df.loc[
            (df['measure_id'] == measure) &
            (df['counterfactual'] == 'no_infection') &
            (df['abx_set'] == 'all') & 
            (df['abx_class'].isin(['all_susceptible', 'all_resistant'])),
            'point_estimate'
        ].sum()
        assert np.isclose(before_check, after_check)
        if 'detailed_location' in df.columns:
            df.loc[df['detailed_location'].notnull(), 'location_id'] = df['detailed_location']

        return df

    def aggregate_2021_locs_to_2023(self, df):
        df = add_location_metadata(df, ['level', 'iso3', 'parent_id'], **self.cache_kwargs)
        df.loc[(df['level'] == 6) & (df['iso3'] == 'GBR'), 'location_id'] = df['parent_id']
        df.loc[(df['level'] == 4) & (df['iso3'] == 'SWE'), 'location_id'] = df['parent_id']
        df = df.groupby([col for col in df.columns if col not in self.draw_cols + ['point_estimate']], as_index=False)\
            [self.draw_cols + ['point_estimate']].sum()
        df.drop(columns=['level', 'iso3', 'parent_id'], inplace=True)
        return df

    def run_calculator(self):
        print_log_message(
            f"Reading in {self.burden_type} counts for pathogen {self.pathogen}, "
            f"year {self.year_id}, cause {self.cause_id}, infectious_syndrome "
            f"{self.infectious_syndrome}"
        )
        df = AmrResult(
            process='split_pathogen',
            burden_type=self.burden_type,
            year_id=self.year_id,
            cause_id=self.cause_id,
            infectious_syndrome=self.infectious_syndrome,
            pathogen=self.pathogen
        ).read_results()
        df = self.aggregate_2021_locs_to_2023(df)

        print_log_message("Converting to year-based metrics")
        if self.burden_type == 'fatal':
            df = df.assign(measure_id=1).append(
                self.calculate_ylls(df).assign(measure_id=4)
            )
        elif self.burden_type == 'nonfatal':
            ypc = AmrResult(
                process='calculate_ylds_per_case', burden_type=self.burden_type,
                year_id=self.year_id, infectious_syndrome=self.infectious_syndrome
            ).read_results()
            if self.infectious_syndrome == 'L2_blood_stream_infection':
                df = add_age_metadata(
                    df, 'age_group_days_end', **self.cache_kwargs
                )
                report_if_merge_fail(df, 'age_group_days_end', 'age_group_id')
                neonatal = df.query("age_group_days_end <= 27")
                non_neonatal = df.query("age_group_days_end > 27")
                if len(neonatal) > 0:
                    neonatal = self.calculate_ylds(
                        neonatal, ypc.loc[ypc.cause_id == 383].copy()
                    )
                if len(non_neonatal) > 0:
                    non_neonatal = self.calculate_ylds(
                        non_neonatal, ypc.loc[ypc.cause_id == 368].copy()
                    )
                if len(df) == 0:
                    AssertionError, 'Empty Dataframe!'
                df = df.assign(measure_id=6).append(
                    pd.concat([neonatal, non_neonatal], sort=False)
                    .assign(measure_id=3)
                )
            else: 
                if self.infectious_syndrome == 'L2_skin_infection':
                    if self.cause_id in [656, 657, 665]:
                        ypc = ypc.query(f"cause_id == {self.cause_id}")
                    else:
                        ypc = ypc.query(f"cause_id == 980")
                elif self.infectious_syndrome == 'L2_typhoid_paratyphoid_ints':
                    ypc = ypc.loc[ypc['cause_id'] == self.cause_id]
                elif self.infectious_syndrome == 'L2_sexually_transmitted_infection':
                    if self.cause_id in [394, 395, 396]:
                        ypc = ypc.loc[ypc['cause_id'] == self.cause_id, ]
                    else:
                        ypc = ypc.loc[ypc['cause_id'] == 393, ]
                elif self.infectious_syndrome == 'L2_hepatitis':
                    if self.cause_id in [401, 402, 403, 404]:
                        ypc = ypc.loc[ypc['cause_id'] == self.cause_id, ]
                    else:
                        ypc = ypc.loc[ypc['cause_id'] == 400, ]
                df = df.assign(measure_id=6).append(
                    self.calculate_ylds(df, ypc).assign(measure_id=3)
                )

        print_log_message("Calculating AMR burden for 2 counterfactuals")
        abxs = self.csbg.loc[(self.csbg['pathogen'] == self.pathogen),
                             'abx_class'].unique().tolist()
        if ('none_tested' not in abxs):
            props = AmrResult(
                process='calculate_amr_props',
                burden_type=self.burden_type,
                year_id=self.year_id,
                pathogen=self.pathogen
            ).read_results()
            if burden_type == 'fatal':
                if self.infectious_syndrome in \
                    ['L2_lower_respiratory_infection', 'L2_blood_stream_infection', 'L2_urinary_tract_infection', 'L2_tb']:
                    props = props.loc[props['infectious_syndrome'] == self.infectious_syndrome, ]
                else:
                    props = props.loc[props['infectious_syndrome'] == 'other_syndromes', ]
            else:
                props['infectious_syndrome'] = 'all'
            props = props.drop(columns=['infectious_syndrome'])
            print_log_message("Got the amr props, multiplying them to split pathogen results")
            df = self.multiply_prop(df, props)
        else:
            print_log_message("No AMR burden for this pathogen, report all pathogen deaths as is")
            df['abx_set'] = 'all'
            df['abx_class'] = '(none_tested)'
            df['counterfactual'] = 'no_infection'

        df = df.loc[~((df['counterfactual'] == 'susceptible_infection') & (df['measure_id'] == 6)), ]

        df = df[
            self.dem_cols + [
                'cause_id', 'infectious_syndrome',
                'pathogen', 'abx_set', 'abx_class', 'hosp',
                'measure_id', 'counterfactual'
            ] + self.draw_cols + ['point_estimate']
        ]
        assert df.notnull().values.all()

        print_log_message("Saving output")
        AmrResult(
            'calculate_amr_burden',
            self.burden_type,
            self.year_id,
            self.cause_id,
            self.infectious_syndrome,
            self.pathogen
        ).write_results(df)
        return df


if __name__ == '__main__':
    burden_type = str(sys.argv[1])
    task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        base_dir = "FILEPATH"
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        cause_id = int(task_row['cause_id'])
        infectious_syndrome = str(task_row['infectious_syndrome'])
        pathogen = str(task_row['pathogen'])
    else:
        year_id = int(sys.argv[2])
        cause_id = int(sys.argv[3])
        infectious_syndrome = str(sys.argv[4])
        pathogen = str(sys.argv[5])
    calc = AmrBurdenCalculator(
        burden_type=burden_type, year_id=year_id, cause_id=cause_id,
        infectious_syndrome=infectious_syndrome, pathogen=pathogen
    )
    calc.run_calculator()
