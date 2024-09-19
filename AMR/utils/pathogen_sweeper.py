import pandas as pd
import numpy as np
from cod_prep.utils import print_log_message, report_if_merge_fail
from mcod_prep.utils.causes import get_all_related_syndromes

class PathogenSweeper():

    def __init__(self, sweep_type, drop_unuseables=True):
        self.pagm = pd.read_csv("FILEPATH")
        self.ish = pd.read_csv("FILEPATH", encoding='latin1')
        self.sweep_type=sweep_type
        self.drop_unuseables = drop_unuseables
        assert self.sweep_type in ['map_aggregates', 'bridge_map']

    def merge_map_and_validate(self, df):
        df = df.merge(self.pagm, on='pathogen', how='left', validate='many_to_one')
        df = self.deal_with_stragglers(df)

        non_merges = df.loc[
            (df['final_list'].isna()) &
            ~(df['pathogen'].isin(self.pagm['final_list'].unique().tolist())),
        ]
        if len(non_merges) > 0:
            print_log_message('Failed pathogen merges to aggregated form..')
            print(non_merges[['pathogen']].drop_duplicates())
            raise AssertionError

        df.loc[df['final_list'].isna(), 'final_list'] = df['pathogen']

        return df

    def deal_with_stragglers(self, df):
        stragglers = {
            'measle_virus': 'measles',
            'dengue_virus': 'dengue',
            'staphylococcus_others': 'staphylococcus_others_coagulase_negative',
        }
        df.loc[df['pathogen'].isin(stragglers.keys()), 'final_list'] = df['pathogen'].map(stragglers)

        return df

    def remerge_additional_cols(self, df):
        existing_cols = list(set(['aggregation_method', 'package_number', 'pathogen_type']).intersection(set(df.columns)))
        df.drop(columns=existing_cols, inplace=True)
        fl = self.pagm[['final_list', 'aggregation_method', 'package_number', 'pathogen_type']].drop_duplicates()
        df = df.merge(fl, on='final_list', how='left', validate='many_to_one')
        report_if_merge_fail(df, 'aggregation_method', 'final_list')

        return df

    def bridge_map(self, df):
        df = self.remerge_additional_cols(df)
        df.loc[df['aggregation_method'] == 'after RDP move to GRAM_NEGATIVE_OTHER', 'final_list'] = 'gram_negative_others'
        df.loc[df['aggregation_method'] == 'after RDP move to GRAM_POSITIVE_OTHER', 'final_list'] = 'gram_positive_others'

        diarrhea_syndromes = get_all_related_syndromes("L2_diarrhea", self.ish)
        etepec = ['enterotoxigenic_escherichia_coli', 'enteropathogenic_escherichia_coli']
        df.loc[~(df['infectious_syndrome'].isin(diarrhea_syndromes)) & (df['pathogen'].isin(etepec)), 'final_list'] = 'escherichia_coli'

        df.loc[df['pathogen'] == 'salmonella_typhi', 'infectious_syndrome'] = 'L3_typhoid_fever'
        df.loc[df['pathogen'] == 'salmonella_paratyphi', 'infectious_syndrome'] = 'L3_paratyphoid_fever'
        df.loc[df['pathogen'] == 'salmonella_ints_non_typhi/paratyphi', 'infectious_syndrome'] = 'L3_ints'

        return df

    def get_computed_dataframe(self, df):
        if 'pathogen_type' in df.columns:
            df.drop(columns=['pathogen_type'], inplace=True)

        if self.sweep_type == 'map_aggregates':
            start_cases = df['cases'].sum()
            if 'deaths' in df.columns:
                start_deaths = df['deaths'].sum()

            df = self.merge_map_and_validate(df)
            df = self.remerge_additional_cols(df)
            
            assert np.isclose(df['cases'].sum(), start_cases)
            if 'deaths' in df.columns:
                assert np.isclose(df['deaths'].sum(), start_deaths)

        elif self.sweep_type == 'bridge_map':
            df = self.bridge_map(df)

        if self.drop_unuseables:
            if 'source' in df.columns:
                df = df.loc[
                    ((df['source'] == 'SOURCE') & (df['pathogen'] == 'other')) |
                    (df['aggregation_method'] != 'drop not useable'),
                ]

        return df
