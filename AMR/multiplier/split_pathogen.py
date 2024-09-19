
import pandas as pd
import sys
import os
import numpy as np
from cod_prep.claude.configurator import Configurator
from cod_prep.utils import print_log_message, report_if_merge_fail, report_duplicates
from cod_prep.downloaders import (
    getcache_age_aggregate_to_detail_map,
    get_current_location_hierarchy,
    prep_child_to_available_parent_loc,
    get_current_cause_hierarchy,
    add_age_metadata,
    add_location_metadata
)
from amr_prep.utils.amr_io import AmrResult
from amr_prep.utils.misc import get_prepped_csbg_universe

CONF = Configurator()
CSV_ID = CONF.get_id('computation_cause_set_version')
CS_ID = CONF.get_id('computation_cause_set')


class PathogenSplitter(object):

    def __init__(self, burden_type, year_id, cause_id, infectious_syndrome):
        self.dem_cols = [
            'location_id', 'sex_id', 'year_id', 'age_group_id', 'hosp'
        ]
        self.cache_kwargs = {'force_rerun': False, 'block_rerun': True, 'cache_results': False}
        self.conf = Configurator()
        self.draw_cols = ["draw_" + str(x) for x in range(0, self.conf.get_id('draw_num'))]
        self.burden_type = burden_type
        self.year_id = year_id
        self.cause_id = cause_id
        self.infectious_syndrome = infectious_syndrome

    def convert_to_incidence(self, df, mir):
        """Calculate incidence cases from deaths for nonfatal burden"""
        print_log_message("Prepping for merge")
        merge_cols = [
            'location_id', 'year_id', 'age_group_id', 'sex_id',
            'infectious_syndrome'
        ]
        if 'hosp' in mir:
            if mir['hosp'].isna().values.all():
                mir.drop(columns=['hosp'], inplace=True)
            else:
                merge_cols += ['hosp']

        if not (set(df.location_id) > set(mir.location_id)):
            df = add_location_metadata(df, 'iso3')
            lh = get_current_location_hierarchy(**self.cache_kwargs)
            iso3_map = lh.loc[lh['level'] == 3, ['iso3', 'location_id']].drop_duplicates()
            iso3_map.rename(columns={'location_id': 'country_loc_id'}, inplace=True)
            df = df.merge(iso3_map, on=['iso3'], how='left', validate='many_to_one')
            assert df['country_loc_id'].notnull().values.all()
            mir.rename(columns={'location_id': 'country_loc_id'}, inplace=True)
            merge_cols.remove('location_id')
            merge_cols.append('country_loc_id')
        if not (set(df.age_group_id) <= set(mir.age_group_id)):
            age_map = getcache_age_aggregate_to_detail_map(**self.cache_kwargs)
            age_map = age_map.loc[age_map.agg_age_group_id.isin(mir.age_group_id.unique())]
            df = df.merge(
                age_map, how='left', on='age_group_id', validate='many_to_one'
            )
            report_if_merge_fail(df, 'agg_age_group_id', 'age_group_id')
            merge_cols.remove('age_group_id')
            merge_cols.append('agg_age_group_id')
            mir = mir.rename(columns={'age_group_id': 'agg_age_group_id'})
        if not (set(df.sex_id) <= set(mir.sex_id)):
            df['parent_sex_id'] = 3
            merge_cols.remove('sex_id')
            merge_cols.append('parent_sex_id')
            mir = mir.rename(columns={'sex_id': 'parent_sex_id'})

        print_log_message("Merging")
        df = df.merge(
            mir, how='outer', on=merge_cols, validate='many_to_one',
            indicator=True, suffixes=('', '_ratio')
        )
        assert (df._merge != 'left_only').all()
        df = df.loc[df._merge != 'right_only']

        print_log_message("Calculating incident cases")
        df[self.draw_cols + ['point_estimate']] = pd.DataFrame(
            df[self.draw_cols + ['point_estimate']].to_numpy() * df[['ratio']].to_numpy(),
            index=df.index
        )

        if self.infectious_syndrome in \
            ['L2_hepatitis', 'L2_tb', 'L2_sexually_transmitted_infection', 'L2_typhoid_paratyphoid_ints']:
            assert set(df.loc[df.hosp.notnull(), 'hosp']) == {'all'}
            df['hosp'] = df['hosp'].fillna('all')
            if ~(df.incidence_only == False).values.all():
                df.loc[df.incidence_only, self.draw_cols] = df.filter(
                    regex='draw_\\d{1,3}_incidence'
                ).rename(columns=lambda x: x.replace("_incidence", ""))
                df.loc[df.incidence_only, 'point_estimate'] = df['pe_incd']
            
        df = df[
            self.dem_cols + self.draw_cols + ['point_estimate', 'infectious_syndrome', 'cause_id']
        ]
        assert df.notnull().values.all()
        assert (df != np.inf).values.all()
        return df

    @staticmethod
    def read_pathogen_models(conf, infectious_syndrome, year_id, **cache_kwargs):
        model_versions = pd.read_csv("FILEPATH")
        model_versions = model_versions.groupby('infectious_syndrome').apply(
            lambda d: dict(zip(d['model_version'], d['query_string']))).to_dict()

        combined_syndromes = pd.read_csv("FILEPATH")
        if infectious_syndrome in combined_syndromes['component_syndromes'].unique().tolist():
            pull_syndrome = combined_syndromes.loc[
                combined_syndromes['component_syndromes'] == infectious_syndrome, 
                'path_dist_agg_syndrome'
            ].item()
        elif infectious_syndrome in ['L2_endocarditis', 'L1.5_myocarditis_pericarditis_carditis']:
            pull_syndrome = 'L1_cardiovascular_infection'
        else:
            pull_syndrome = infectious_syndrome

        bug = pd.DataFrame()
        for mv, query_string in model_versions[pull_syndrome].items():
            model_df = pd.read_csv("FILEPATH")
            if type(query_string) == str:
                model_df = add_age_metadata(
                    model_df, ['age_group_days_start', 'age_group_days_end'],
                    **cache_kwargs
                )
                report_if_merge_fail(model_df, 'age_group_days_start', 'age_group_id')
                model_df = model_df.query(query_string)
                model_df = model_df.drop(
                    ['age_group_days_start', 'age_group_days_end'],
                    axis='columns'
                )
            bug = bug.append(model_df)
        bug['infectious_syndrome'] = infectious_syndrome
        return bug

    def split_pathogen(self):
        print_log_message(f"Reading pathogen draws")

        bug = self.read_pathogen_models(
            self.conf, self.infectious_syndrome, self.year_id, **self.cache_kwargs)

        print_log_message("Multiplying by pathogen distribution")
        if self.burden_type == 'fatal':
            bug = bug.loc[bug.measure_id == 1]
        elif self.burden_type == 'nonfatal':
            bug = bug.loc[bug.measure_id == 6]

        if (self.infectious_syndrome == 'L2_diarrhea'):
            pathogen_merge_cols = [
                'infectious_syndrome', 'location_id',
                'age_group_id', 'sex_id', 'year_id'
            ]
            diarrhea_rename = {
                "vibrio_cholerae_spp": "vibrio_cholerae",
                "non_typhoidal_salmonellae": "salmonella_ints_non_typhi",
                "campylobacter": "campylobacter_spp",
                "amebiasis": "entamoeba_histolytica",
                "cryptosporidiosis": "cryptosporidium_spp",
                "adenoviruses": "adenovirus"
            }
            bug['pathogen_rename'] = bug['pathogen'].map(diarrhea_rename)
            bug.loc[bug['pathogen_rename'].notnull(), 'pathogen'] = bug['pathogen_rename']
            bug.drop(columns=['pathogen_rename'], inplace=True)
        else:
            age_map = getcache_age_aggregate_to_detail_map(**self.cache_kwargs)
            age_map = age_map.loc[age_map.agg_age_group_id.isin(bug.age_group_id.unique())]
            self.df = self.df.merge(age_map, how='left', on='age_group_id', validate='many_to_one')
            report_if_merge_fail(self.df, 'agg_age_group_id', 'age_group_id')
            if self.infectious_syndrome == 'L2_urinary_tract_infection':
                self.df['parent_sex_id'] = self.df['sex_id']
                bug['sex_id'] = bug['sex_id'].map({0:1, 1:2})
            else:
                self.df['parent_sex_id'] = 3
            lh = get_current_location_hierarchy(
                location_set_version_id=self.conf.get_id("location_set_version"),
                **self.cache_kwargs)
            loc_map = prep_child_to_available_parent_loc(
                lh, bug.location_id.unique().tolist(), min_level=3)
            self.df = self.df.merge(loc_map, how='left', on='location_id', validate='many_to_one')
            report_if_merge_fail(self.df, 'parent_location_id', 'location_id')
            bug = bug.rename(
                columns={
                    'age_group_id': 'agg_age_group_id',
                    'location_id': 'parent_location_id',
                    'sex_id': 'parent_sex_id'
                }
            )
            pathogen_merge_cols = [
                'agg_age_group_id', 'parent_sex_id', 'parent_location_id',
                'year_id', 'infectious_syndrome', 'hosp'
            ]
        assert np.allclose(bug.groupby(pathogen_merge_cols)[self.draw_cols].sum(), 1)
        bug['point_estimate'] = bug[self.draw_cols].mean(1)
        report_duplicates(bug, pathogen_merge_cols + ['pathogen'])

        self.df = self.df.merge(
            bug, how='left',
            on=pathogen_merge_cols,
            suffixes=('', '_path')
        )
        report_if_merge_fail(self.df, self.draw_cols[0] + '_path', pathogen_merge_cols)
        self.df[self.draw_cols] = pd.DataFrame(
            self.df[self.draw_cols].to_numpy() * self.df.filter(
                regex='draw_\\d{1,3}_path').to_numpy(),
            index=self.df.index
        )
        self.df['point_estimate'] = self.df['point_estimate'] * self.df['point_estimate_path']
        self.df = self.df.drop([d + '_path' for d in self.draw_cols] + ['point_estimate_path'], axis='columns')

    def split_others(self):
        print_log_message('Splitting "other" pathogen into specific remainders')
        prop_other_dir = CONF.get_resource('prop_for_others')
        if self.infectious_syndrome == 'L2_endocarditis':
            po = pd.read_csv("FILEPATH")
        else:
            po = pd.read_csv("FILEPATH")
        df_other = self.df.loc[self.df['pathogen'] == 'other', ]
        self.df = self.df.loc[self.df['pathogen'] != 'other', ]
        df_other.rename(columns={'pathogen': 'presplit_pathogen'}, inplace=True)
        po['presplit_pathogen'] = 'other'
        before_sum = df_other['point_estimate'].sum()
        df_other = df_other.merge(po, on='presplit_pathogen', how='outer')
        df_other.reset_index(inplace=True)
        df_other[self.draw_cols + ['point_estimate']] = pd.DataFrame(
            df_other[self.draw_cols + ['point_estimate']].to_numpy() * df_other[['prop']].to_numpy(),
            index=df_other.index
        )
        df_other.drop(columns=['prop', 'presplit_pathogen'], inplace=True)
        assert np.isclose(df_other['point_estimate'].sum(), before_sum)
        self.df = pd.concat([self.df, df_other])

    def run_split(self):
        print_log_message(
            "Working on " + str(self.year_id) + ' ' + str(self.cause_id) + ' ' + self.infectious_syndrome
        )
        print_log_message(f"Reading syndrome draws")
        self.df = AmrResult(
            process='split_sepsis_syndrome',
            burden_type='fatal',
            year_id=self.year_id,
            cause_id=self.cause_id,
            infectious_syndrome=self.infectious_syndrome
        ).read_results()

        if self.burden_type == 'nonfatal':
            print_log_message("Reading MI ratios")
            mir = AmrResult(
                process='calculate_mi_ratios',
                burden_type=self.burden_type,
                year_id=self.year_id,
                infectious_syndrome=self.infectious_syndrome
            ).read_results()
            no_cfr_syndromes = [
                'L2_endocarditis', 'L1.5_myocarditis_pericarditis_carditis', 'L2_tb',
                'L2_upper_respiratory_infection', 'L2_hepatitis', 'L2_diarrhea',
                'L2_sexually_transmitted_infection', 'L2_typhoid_paratyphoid_ints',
                'L2_other_parasitic_infection', 'L2_unspecified_site_infection', 'L2_OAMOIS'
            ]
            if self.infectious_syndrome in no_cfr_syndromes:
                if self.cause_id in [319, 320, 959, 401, 402, 403, 404, 394, 395, 396]:
                    mir = mir.loc[mir['cause_id'] == self.cause_id, ]
                elif (self.infectious_syndrome == 'L2_hepatitis') & \
                    (self.cause_id not in [401, 402, 403, 404]):
                    mir = mir.loc[mir['cause_id'] == 400, ]
                elif (self.infectious_syndrome == 'L2_sexually_transmitted_infection') & \
                    (self.cause_id == 399):
                    mir = mir.loc[mir['cause_id'] == 393, ]
                mir.rename(columns={'mi_ratio': 'ratio'}, inplace=True)
                self.df = self.convert_to_incidence(self.df, mir)
            else:
                if 'cause_id' in mir.columns:
                    mir = mir.loc[mir.cause_id.isnull(), ]
                mir.rename(columns={'1/cfr': 'ratio'}, inplace=True)
                self.df = self.convert_to_incidence(self.df, mir)

        csbg_universe = get_prepped_csbg_universe(
            burden_type=self.burden_type, add_parents=False
        )
        csb = csbg_universe.loc[
            (csbg_universe['cause_id'] == self.cause_id) &
            (csbg_universe['infectious_syndrome'] == self.infectious_syndrome),
            ['cause_id', 'infectious_syndrome', 'pathogen']
        ].drop_duplicates()
        no_path_syns = csb.loc[csb['pathogen'].isna(), 'infectious_syndrome'].unique().tolist()
        if (len(csb) == 1) & (csb['pathogen'].notnull().any()):
            print_log_message("This syndrome/cause is 1 pathogen")
            pathogen = csb['pathogen'].unique()[0]
            self.df['pathogen'] = pathogen
        elif self.infectious_syndrome == 'L2_tb':
            self.df['pathogen'] = 'mycobacterium_tuberculosis'
        elif self.infectious_syndrome in no_path_syns:
            self.df['pathogen'] = '(none_estimated)'
        else:
            self.split_pathogen()

        syndromes_split_other = [
            'L2_meningitis',
            'L1_bone_joint_infection',
            'L2_endocarditis',
            'L1_peritoneal_and_intra_abdomen_infection',
            'L2_blood_stream_infection',
            'L2_lower_respiratory_infection',
            'L2_skin_infection',
            'L2_urinary_tract_infection'
        ]
        if self.infectious_syndrome in syndromes_split_other:
            self.split_others()

        self.df = self.df[
            self.dem_cols + ['cause_id', 'infectious_syndrome', 'pathogen']
            + self.draw_cols + ['point_estimate']
        ]
        assert self.df.notnull().values.all()

    def save_output(self):
        print_log_message("Saving outputs")
        for pathogen, pathogen_df in self.df.groupby('pathogen'):
            print_log_message(f"Working on pathogen {pathogen}")
            AmrResult(
                process='split_pathogen',
                burden_type=self.burden_type,
                year_id=self.year_id,
                cause_id=self.cause_id,
                infectious_syndrome=self.infectious_syndrome,
                pathogen=pathogen
            ).write_results(pathogen_df)


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
    else:
        year_id = int(sys.argv[2])
        cause_id = int(sys.argv[3])
        infectious_syndrome = str(sys.argv[4])

    splitter = PathogenSplitter(
        burden_type=burden_type,
        year_id=year_id, cause_id=cause_id,
        infectious_syndrome=infectious_syndrome
    )
    splitter.run_split()
    splitter.save_output()
