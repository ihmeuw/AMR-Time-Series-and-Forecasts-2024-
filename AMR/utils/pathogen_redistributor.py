import pandas as pd
import numpy as np
from cod_prep.downloaders.locations import *
from cod_prep.downloaders import *
from cod_prep.utils import print_log_message, report_if_merge_fail
from amr_prep.utils.pathogen_sweeper import PathogenSweeper
from cod_prep.claude.configurator import Configurator


class PathogenRedistributer():
    CONF = Configurator()

    def __init__(
        self,
        drop_remainder=False,
        bridge_map=True,
        revert_col_names=True,
        verbose=False,
        drop_unuseables=False,
        data_type=None,
        add_cutoff=None
    ):
        self.props = pd.read_csv("FILEPATH")
        self.pagm = pd.read_csv("FILEPATH")
        self.ish = pd.read_csv("FILEPATH", encoding='latin1')
        self.year_bins = [(1980, 2004), (2005, 2023)]
        self.agg_age_groups_years_name = {
            (0, 0.07671233): 'Neonatal',
            (0, 1): '<1 year',
            (0, 5): '<5 years',
            (0.07671233, 5): 'Post Neonatal to 5',
            (5, 125): '5 plus'
        }
        self.target_genuses = \
            self.pagm.loc[self.pagm['aggregation_method'] == 'redistributable','final_list'].unique().tolist()
        self.agg_columns_list = [
            'WB Income Level',
            'year_bin',
            'aggregated_age_group',
            'L1_infectious_syndrome'
        ]
        self.collapse_order_values = {
            'WB Income Level': 'All Incomes',
            'year_bin': str(self.year_bins[0][0]) + '-' + str(self.year_bins[1][1]),
            'aggregated_age_group': 'All Ages',
        }
        self.minimum_pathogen_num_prop = 2
        self.drop_remainder = drop_remainder
        self.bridge_map = bridge_map
        self.revert_col_names = revert_col_names
        self.verbose=verbose
        self.drop_unuseables=drop_unuseables
        self.data_type = data_type
        self.add_cutoff = add_cutoff
        self.props.loc[self.props['year_bin'] == '1980-2022', 'year_bin'] = '1980-2023'
        self.props.loc[self.props['year_bin'] == '2005-2022', 'year_bin'] = '2005-2023'
        self.add_cns_bone_non_pneumo_props()
        self.add_lri_non_pneumo_groupB_props()

    def add_cns_bone_non_pneumo_props(self):
        snp = self.props.loc[
            (self.props['infectious_syndrome'].isin(['L1_cns_infection', 'L1_bone_joint_infection'])) & 
            (self.props['genus'] == 'streptococcus') &
            (self.props['pathogen'] != 'streptococcus_pneumoniae'),
        ]
        snp['total_cases_genus'] = snp.groupby(
            ['WB Income Level', 'year_bin', 'age_group_name', 'infectious_syndrome', 'genus'], 
            as_index=False
        )['species_cases'].transform(sum)

        snp['genus'] = 'streptococcus_non_pneumoniae'
        snp['prop_cases'] = snp['species_cases'] / snp['total_cases_genus']

        snp['# species in prop'] = snp.groupby(
            ['WB Income Level', 'infectious_syndrome', 'genus', 'age_group_name', 'year_bin'],
            as_index=False
        )['pathogen'].transform('nunique')

        snp['prop_id'] = snp.groupby(
            ['WB Income Level', 'infectious_syndrome', 'genus', 'age_group_name', 'year_bin']
        ).ngroup().add(1)

        snp['prop_sum'] = snp.groupby('prop_id')['prop_cases'].transform('sum')
        assert np.isclose(snp['prop_sum'], 1).all()
        snp['prop_id'] = snp['prop_id'] + .5
        self.props = pd.concat([self.props, snp])

    def add_lri_non_pneumo_groupB_props(self):
        snp = self.props.loc[
            (self.props['infectious_syndrome'] == 'L1_respiratory_infection') & 
            (self.props['genus'] == 'streptococcus') &
            ~(self.props['pathogen'].isin(['streptococcus_pneumoniae', 'streptococcus_group_b'])),
        ]
        snp['total_cases_genus'] = snp.groupby(
            ['WB Income Level', 'year_bin', 'age_group_name', 'infectious_syndrome', 'genus'], 
            as_index=False
        )['species_cases'].transform(sum)

        snp['genus'] = 'strep_non_pneumo_group_b'
        snp['prop_cases'] = snp['species_cases'] / snp['total_cases_genus']

        snp['# species in prop'] = snp.groupby(
            ['WB Income Level', 'infectious_syndrome', 'genus', 'age_group_name', 'year_bin'],
            as_index=False
        )['pathogen'].transform('nunique')

        snp['prop_id'] = snp.groupby(
            ['WB Income Level', 'infectious_syndrome', 'genus', 'age_group_name', 'year_bin']
        ).ngroup().add(1)

        snp['prop_sum'] = snp.groupby('prop_id')['prop_cases'].transform('sum')
        assert np.isclose(snp['prop_sum'], 1).all()
        snp['prop_id'] = snp['prop_id'] + .75
        self.props = pd.concat([self.props, snp])

    def subset_input_data(self, df):
        aggregator = PathogenSweeper(sweep_type='map_aggregates', drop_unuseables=self.drop_unuseables)
        df = aggregator.get_computed_dataframe(df)

        non_rd = df.loc[
            ~(df['final_list'].isin(self.target_genuses)) |
            (df['year_id'] < 1980),
        ]
        df = df.loc[
            (df['final_list'].isin(self.target_genuses)) &
            (df['year_id'] >= 1980),
        ]

        return df, non_rd

    def add_merge_cols(self, df):
        if self.verbose:
            print_log_message("Merging on aggregated values for location, years, age, and syndrome...")
        wblh = get_current_location_hierarchy(
            release_id=CONF.get_id('release'),
            location_set_id=26,
            location_set_version_id=1115
        )
        wblh.rename(columns={'region_name': 'WB Income Level'}, inplace=True)
        lh = get_current_location_hierarchy()
        lh = lh.merge(wblh[['location_id', 'WB Income Level']], on='location_id', how='left')                    
        lh = lh.loc[lh['WB Income Level'].notnull(), ['iso3', 'WB Income Level']].drop_duplicates()
        df = add_location_metadata(df, 'iso3')
        df = df.merge(lh[['iso3', 'WB Income Level']], on='iso3', how='left', validate='many_to_one')
        df.loc[df['WB Income Level'] != 'World Bank High Income', 'WB Income Level'] = 'World Bank Middle and Low Income'

        for yb in self.year_bins:  
            df.loc[df['year_id'].between(yb[0], yb[1]), 'year_bin'] = str(yb[0]) + '-' + str(yb[1])

        df = add_age_metadata(df, ['age_group_years_start', 'age_group_years_end'])
        for agyb, agyn in self.agg_age_groups_years_name.items():
            df.loc[
                (df['age_group_years_start'] >= agyb[0]) &
                (df['age_group_years_end'] <= agyb[1]),
            'aggregated_age_group'] = agyn
        df.loc[df['aggregated_age_group'].isna(), 'aggregated_age_group'] = 'All Ages'

        is_ptp = self.ish[['infectious_syndrome', 'path_to_top_parent']]
        is_ptp[['0', '1', '2', '3', '4']] = is_ptp['path_to_top_parent'].str.split(',', expand=True)
        df = df.merge(is_ptp, how='left', on='infectious_syndrome', validate='many_to_one')
        df.drop(columns=['0', '2', '3', '4'], inplace=True)
        df.rename(columns={'1': 'L1_infectious_syndrome'}, inplace=True)

        assert df[self.agg_columns_list].notnull().values.all()

        return df

    def merge_props(self, df):
        df['genus'] = df['final_list'].str.replace('_unspecified', '')
        df = df.merge(
            self.props,
            on=self.agg_columns_list + ['genus'],
            how='outer',
            indicator=True,
        )
        df_success = df.loc[df['_merge'] == 'both', ]
        df_fail = df.loc[df['_merge'] == 'left_only', ]
        df_others_fail = df_success.loc[
            (df_success['# species in prop'] == 1) &
            (df_success['genus_species'].str.contains('_others')),
        ]
        df_success = df_success.loc[
            ~((df_success['# species in prop'] == 1) &
              (df_success['genus_species'].str.contains('_others'))),
        ]
        df_fail = pd.concat([df_fail, df_others_fail])
        prop_value_cols = list(set(self.props.columns) - set(self.agg_columns_list + ['genus'])) + ['_merge']
        df_fail.drop(columns=prop_value_cols, inplace=True)

        return df_success, df_fail

    def redistribute_pathogens(self, df):
        if self.verbose:
            print_log_message('Merge and apply proportions...')

        prop_col_rename_dict = dict(
            zip(
                ['WB Income Level', 'year_bin', 'age_group_name', 'infectious_syndrome', 'pathogen'],
                self.agg_columns_list + ['genus_species']
            )
        )
        self.props.rename(columns=prop_col_rename_dict, inplace=True)
        if self.add_cutoff is not None:
            self.props = self.props.loc[self.props['total_cases_genus'] >= self.add_cutoff, ] 

        self.remainder_cases = 0
        self.remainder_deaths = 0
        df_success1, df_fail = self.merge_props(df)
        dfs = []
        for col, val in self.collapse_order_values.items():
            df_fail[col] = val
            df_success2, df_fail = self.merge_props(df_fail)
            dfs.append(df_success2)
            if len(df_fail) == 0:
                break
            if (col == 'aggregated_age_group') & (len(df_fail) > 0):
                if self.verbose:
                    print_log_message(
                        "We exhausted all our collapsible demographic columns but unmerged raw data rows remain, "
                        "should we sacrifice the level 1 syndrome too? No, for now just drop these data"
                    )
                    print_log_message('These data without proper merge. Total of ' + str(df_fail['cases'].sum()) + ' cases.')
                self.remainder_cases = df_fail['cases'].sum()
                self.remainder_deaths = df_fail['deaths'].sum()
        df_success2 = pd.concat(dfs)

        df = pd.concat([df_success1, df_success2])
        df['cases'] = df['cases'] * df['prop_cases']
        df['deaths'] = df['deaths'] * df['prop_cases']
        df = df.reset_index(drop=True)
        df.loc[df['genus_species'].notnull(), 'final_list'] = df['genus_species']

        df.drop(columns=self.props.columns, inplace=True)

        if self.drop_remainder:
            if self.verbose:
                print_log_message('We will not include the data without a proper merge..')
        else:
            if self.verbose:
                print_log_message('Unsuccesful merge remainders will be kept..')
            df = pd.concat([df, df_fail])

        return df

    def get_computed_dataframe(self, df):

        df, non_rd = self.subset_input_data(df)

        start_cases = (non_rd['cases'].sum() + df['cases'].sum())
        if 'deaths' not in df.columns:
            df['deaths'] = np.nan
            non_rd['deaths'] = np.nan
        start_deaths = (non_rd['deaths'].sum() + df['deaths'].sum())

        df['pathogen_redistributed'] = 1
        non_rd['pathogen_redistributed'] = 0

        if len(df) == 0:
            if self.verbose:
                print_log_message('No rows of pathogens to redistribute')
            if self.revert_col_names:
                non_rd.rename(columns={'pathogen': 'pre_aggregate_pathogen', 'final_list': 'pathogen'}, inplace=True)
            return non_rd

        else:
            df = self.add_merge_cols(df)
            df = self.redistribute_pathogens(df)
            df = pd.concat([df, non_rd])
            df.reset_index(inplace=True, drop=True)

            if self.drop_remainder:
                assert np.isclose(df['cases'].sum(), start_cases - self.remainder_cases)
                assert np.isclose(df['deaths'].sum(), start_deaths - self.remainder_deaths)
            else:
                assert np.isclose(df['cases'].sum(), start_cases)
                assert np.isclose(df['deaths'].sum(), start_deaths)

            if self.drop_remainder:
                assert len(df.loc[(df['final_list'].isin(self.target_genuses)) & (df['year_id'] >= 1980), ]) == 0
            
            if self.verbose:
                print_log_message("Finished redistributing pathogens..")

            if self.bridge_map:
                if self.verbose:
                    print_log_message('Bridge mapping redistribution artifacts..')
                bridge_mapper = PathogenSweeper('bridge_map', drop_unuseables=self.drop_unuseables)
                df = bridge_mapper.get_computed_dataframe(df)

            if self.verbose:
                print_log_message('Dropping any pre 1980 data')
            df = df.loc[df['year_id'] >= 1980, ]

            if self.revert_col_names:
                df.rename(columns={'pathogen': 'pre_aggregate_pathogen', 'final_list': 'pathogen'}, inplace=True)

            if df['deaths'].isnull().values.all():
                df.drop(columns=['deaths'], inplace=True)

            return df
