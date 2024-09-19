import pandas as pd
import numpy as np
from warnings import warn
from cod_prep.claude.claude_io import Configurator
from mcod_prep.utils.causes import get_infsyn_hierarchy, get_all_related_syndromes
from mcod_prep.utils.mcause_io import get_mcause_data
from amr_prep.utils.amr_io import get_amr_data
from cod_prep.utils import (
    print_log_message, report_duplicates, create_square_df,
    report_if_merge_fail
)
from mcod_prep.utils.nids import add_nid_metadata, get_datasets
from cod_prep.claude.relative_rate_split import relative_rate_split
from cod_prep.downloaders import (
    add_age_metadata, get_current_location_hierarchy,
    add_location_metadata, get_cod_ages, pretty_print,
    get_pop, add_population, getcache_age_aggregate_to_detail_map,
    get_country_level_location_id, get_ages,
    prep_age_aggregate_to_detail_map
)
from amr_prep.utils.pathogen_redistributor import PathogenRedistributer


class PathogenFormatter():
    blank_pathogens = ['none', 'unknown']
    unique_cols = [
        'source', 'nid', 'data_type', 'location_id', 'year_id', 'age_group_id', 'sex_id', 'pre_aggregate_pathogen',
        'hosp', 'infectious_syndrome', 'main_diagnosis', 'sample_id', 'pathogen', 'pathogen_type',
        'pathogen_redistributed',
    ]
    flu_rsv_surv = ["SOURCES"]

    def __init__(self, model_type, infectious_syndrome, keep_pathogens,
                 age_weights_use=None, cache_kwargs=None, infsyn_version=None):
        self.model_type = model_type
        self.infsyn_version = infsyn_version
        self.infectious_syndrome = infectious_syndrome
        self.keep_pathogens = keep_pathogens
        self.age_weights_use = age_weights_use or {}
        if not cache_kwargs:
            self.cache_kwargs = {
                'force_rerun': False, 'block_rerun': True, 'cache_results': False
            }
        else:
            self.cache_kwargs = cache_kwargs

        self.conf = Configurator()
        self.validate_inputs()
        self.set_sources()

    def set_sources(self):
        self.amr_metadata = pd.read_csv("FILEPATH")
        self.mcause_metadata = get_datasets(
            is_active=True, data_type_id=[3, 9, 13, 16],
            **self.cache_kwargs
        )

        if self.model_type in ['network', 'flu_rsv_cw']:
            amr_source_query = "active_3a == 1"
            mcause_data_type_query = "data_type_id in (3, 9, 13, 16)"
        elif self.model_type == 'cfr':
            amr_source_query = "active_3b == 1"
            mcause_data_type_query = "data_type_id in (3, 16)"
        self.pull_sources = set(self.amr_metadata.query(
            amr_source_query + " | inform_age_weights == 1"
        ).source)
        self.final_sources = set(
            self.amr_metadata.query(amr_source_query).source
        ).union(set(
            self.mcause_metadata.query(mcause_data_type_query).source
        ))

        self.age_weight_sources = set(self.amr_metadata.query(
            "inform_age_weights == 1"
        ).source).union(self.mcause_metadata.source)

        self.other_node_sources = set(self.amr_metadata.query(
            "inform_other_node == 1"
        ).source)

        self.deaths_only_sources = set(
            self.amr_metadata.query("deaths_only == 1").source
        ).union(
            self.mcause_metadata.query("data_type_id in (9, 13)").source
        )

        self.icu_only_sources = set(
            self.amr_metadata.query("icu_only == 1").source
        )


    def validate_inputs(self):
        assert self.conf.config_type == 'amr'
        assert self.model_type in ['network', 'cfr', 'flu_rsv_cw']

        infsyn = get_infsyn_hierarchy(infsyn_set_version=self.infsyn_version)
        self.infsyn = infsyn 
        combined_syndromes = pd.read_csv("FILEPATH")
        if self.infectious_syndrome not in self.infsyn.infectious_syndrome.unique():
           assert self.infectious_syndrome in combined_syndromes['path_dist_agg_syndrome'].unique()
        if self.infectious_syndrome in combined_syndromes['path_dist_agg_syndrome'].unique().tolist():
            component_syns = combined_syndromes.loc[
                combined_syndromes['path_dist_agg_syndrome'] == self.infectious_syndrome,
                'component_syndromes'
            ].unique()
            self.target_syndromes = []
            for syn in component_syns:
                self.target_syndromes += get_all_related_syndromes(syn, self.infsyn)
        else:
            self.target_syndromes = get_all_related_syndromes(
                self.infectious_syndrome, self.infsyn)
        if self.model_type == 'flu_rsv_cw':
            infsyn_v2 = get_infsyn_hierarchy('v2_ari')
            self.target_syndromes = list(set(self.target_syndromes).union(
                get_all_related_syndromes('acute_respiratory_infectious', infsyn_v2)
            ))

        final_pathogens = pd.read_csv("FILEPATH")
        final_pathogens = final_pathogens.loc[final_pathogens['modelling_side'] != 'GBD', ]
        assert set(self.keep_pathogens) <= set(final_pathogens.pathogen)
        assert set(self.keep_pathogens).isdisjoint({'other'})
        assert set(self.age_weights_use) <= set(final_pathogens.pathogen).union(
            {'all', 'other'})
        assert set(self.age_weights_use.values()) <= set(final_pathogens.pathogen).union(
            {'all', 'other'})
        if self.model_type == 'flu_rsv_cw':
            assert (set(self.keep_pathogens) == {'flu'}) or (
                set(self.keep_pathogens) == {'rsv'})

    def fix_ages(self, df):
        df.loc[df.age_group_id == 325, 'age_group_id'] = 238
        df.loc[df.age_group_id.isin([247, 361, 382]), 'age_group_id'] = 183
        df.loc[df.age_group_id.isin([164, 429]), 'age_group_id'] = 2
        return df

    def drop_data_for_flu_rsv_cw(self, df):
        potential_drops = df.loc[
            ~df.data_type.isin(['mcod', 'hospital', 'linkage'])
            & ~df.source.isin(self.other_node_sources), 'source'].unique()

        flu_rsv_full_denom = ["SOURCES"]
        potential_drops = set(potential_drops) - set(flu_rsv_full_denom)

        problem_sources = df.loc[
            df.source.isin(potential_drops)
            & ~df.source.isin(["SOURCES"])
            & df.pathogen.isin(["SOURCES"]),
            'source'
        ].unique().tolist()
        if len(problem_sources) > 0:
            raise AssertionError(
                f"Sources {problem_sources} have data on flu/RSV but do not have a"
                f" full denominator of LRI or ARI - we need a method of dealing with"
                f" this in order to use the data in the crosswalk"
            )

        df = df.loc[~df.source.isin(potential_drops)]
        return df

    def get_data(self):
        print_log_message(
            f"Reading MCoD, hospital, and linkage data"
        )   
        df_mcause = get_mcause_data(
            'format_map', is_active=True, data_type_id=[3, 9, 13, 16],
            sub_dirs='infectious_syndrome', **self.cache_kwargs
        )
        df_mcause = df_mcause.loc[df_mcause['nid'] != 494702, ]
        df_mcause = add_nid_metadata(
            df_mcause, ['source', 'data_type_id'], **self.cache_kwargs)

        df_mcause['infectious_syndrome'] = df_mcause['infectious_syndrome'].astype(str)
        df_mcause = df_mcause.loc[df_mcause.infectious_syndrome.isin(self.target_syndromes)]

        df_mcause['pathogen_from_cause'] = df_mcause.pathogen_from_cause.str.split(',')
        df_mcause['sample_id'] = ['mcause_' + str(x) for x in range(0, len(df_mcause))]
        df_mcause = df_mcause.explode('pathogen_from_cause')
        df_mcause = df_mcause.rename(columns={'pathogen_from_cause': 'pathogen'})
        df_mcause = df_mcause.rename(columns={
            'admissions': 'cases',
            'cause_infectious_syndrome': 'main_diagnosis'
        })
        if 'cases' not in df_mcause.columns:
            df_mcause['cases'] = df_mcause['deaths']
        else:
            df_mcause['cases'] = df_mcause['cases'].fillna(df_mcause['deaths'])

        redistributer = PathogenRedistributer(drop_remainder=True, drop_unuseables=True, add_cutoff=50)
        df_mcause.loc[
            df_mcause['infectious_syndrome'].isin(
                ['bsi_gram_negative_other', 'lri_gram_negative_other']
            ) &
            (df_mcause['pathogen'] == 'gram_negative_other'),
            'pathogen'
        ] = 'gram_negative_unspecified'
        if self.model_type == 'cfr':
            df_mcause_claims = df_mcause.loc[df_mcause['data_type_id'] == 16, ]
            df_mcause = df_mcause.loc[df_mcause['data_type_id'] != 16]
            df_mcause = redistributer.get_computed_dataframe(df_mcause)
            redistributer_claims = PathogenRedistributer(drop_remainder=True, drop_unuseables=False)
            df_mcause_claims = redistributer_claims.get_computed_dataframe(df_mcause_claims)
            df_mcause = pd.concat([df_mcause, df_mcause_claims])
            df_mcause.reset_index(drop=True, inplace=True)
        else:
            df_mcause = redistributer.get_computed_dataframe(df_mcause)

        print_log_message("Reading AMR data...")
        df_amr = get_amr_data(
            sources=list(self.pull_sources),
            infectious_syndrome=self.target_syndromes,
            add_cutoff=50
        )
        df_amr = df_amr.loc[
            ~((df_amr.source == 'SOURCE')
              & (df_amr.cause_id == 368))
        ]
        if self.infectious_syndrome in get_all_related_syndromes(
            'L1_peritoneal_and_intra_abdomen_infection',
            self.infsyn,
            level=2
        ):
            drop_nids = ["NUMBERS"]
            df_amr = df_amr.loc[~df_amr.nid.isin(drop_nids)]

        df_amr['main_diagnosis'] = 'microbiology'

        df = pd.concat([df_mcause, df_amr], sort=False)
        df['hosp'] = df['hosp'].fillna('unknown')

        nid_meta = pd.read_csv("FILEPATH")
        nid_unmapped_meta = pd.read_csv("FILEPATH")
        data_types = nid_meta.append(nid_unmapped_meta).loc[:, [
            'source', 'data_type']].drop_duplicates().reset_index()
        df = df.merge(data_types, on='source', how='left')
        df['data_type'] = df['data_type'].fillna(
            df['data_type_id'].map({9: 'mcod', 3: 'hospital', 13: 'linkage', 16: 'claims'})
        )
        no_dt_sources = list(df.loc[df['data_type'].isna(), 'source'].unique())
        assert len(no_dt_sources) == 0,\
            f"{no_dt_sources} don't have a data_type in their NID metadata, fix imminently!"
        df = df.drop('data_type_id', axis='columns')

        df = df.reset_index(drop=True)
        tabulated_data = df.sample_id.isnull()
        df.loc[tabulated_data, 'sample_id'] = df['source'] + '_'
        df.loc[tabulated_data, 'sample_id'] += [str(x) for x in range(0, tabulated_data.sum())]
        assert df[self.unique_cols].notnull().values.all()
        inconsistent_hosp = df[['source', 'sample_id', 'hosp']].drop_duplicates()\
            .duplicated(subset=['source', 'sample_id']).sum()
        warn(
            f"{inconsistent_hosp} out of {df.sample_id.nunique()} sample_ids"
            f" show inconsistent hosp"
        )
        inconsistent_counts = df[self.unique_cols + ['deaths', 'cases']].drop_duplicates()\
            .duplicated(subset=self.unique_cols, keep=False).sum()
        df = df.sort_values(by='cases', ascending=False).drop_duplicates(
            subset=self.unique_cols, keep='first'
        )
        warn(
            f"{inconsistent_counts} out of {len(df)} records show inconsistent counts"
        )
        df = df[self.unique_cols + ['deaths', 'cases']]

        df = self.fix_ages(df)

        df = df.loc[df.year_id >= 1975]
        df.loc[df.year_id < 1980, 'year_id'] = 1980

        if self.model_type == 'flu_rsv_cw':
            df = self.drop_data_for_flu_rsv_cw(df)

        df = df.loc[df['year_id'] != 2023, ]

        return df

    def select_pathogens(self, df):
        assert set(PathogenFormatter.blank_pathogens).isdisjoint(self.keep_pathogens)
        if self.model_type != 'flu_rsv_cw':
            df = df.loc[~df.pathogen.isin(PathogenFormatter.blank_pathogens)]
        else:
            df = df.loc[~(
                df.pathogen.isin(PathogenFormatter.blank_pathogens)
                & df.source.isin(self.other_node_sources)
            )]


        parents_to_children = {
            "escherichia_coli": r"(enteropathogenic|enterotoxigenic)_escherichia_coli"
        }
        for parent, children in parents_to_children.items():
            df.loc[df.pathogen.str.contains(children), "pathogen"] = parent
        df = df.drop_duplicates(subset=self.unique_cols)

        if self.infectious_syndrome in get_all_related_syndromes(
            'L1_urogenital_infection',
            self.infsyn,
            level=2
        ):
            df = df.loc[~df.pathogen.isin(['chlamydia_spp', 'neisseria_gonorrheae'])]
        if self.infectious_syndrome in get_all_related_syndromes(
            'L1_respiratory_infection',
            self.infsyn,
            level=2
        ):
            df = df.loc[~df.pathogen.str.contains('malaria'), ]
        if self.infectious_syndrome in get_all_related_syndromes(
            'L1_skin_oral_eyes_infection',
            self.infsyn,
            level=2
        ):
            df = df.loc[df.pathogen_type != 'fungi']
        cons = [
            'staphylococcus_epidermidis_coagulase_negative',
            'staphylococcus_haemolyticus_coagulase_negative',
            'staphylococcus_saprophyticus_coagulase_negative',
            'staphylococcus_hominis_coagulase_negative',
            'staphylococcus_others_coagulase_negative',
        ]
        if self.infectious_syndrome not in ['L2_meningitis', 'L2_urinary_tract_infection']:
            df = df.loc[~df['pathogen'].isin(cons), ]
        elif self.infectious_syndrome == 'L2_meningitis':
            df.loc[df['pathogen'].isin(cons), 'pathogen'] = 'coagulase_negative_staphylococcus'
        else:
            cons.remove('staphylococcus_saprophyticus_coagulase_negative')
            df = df.loc[~df['pathogen'].isin(cons), ]
            df.loc[
                df['pathogen'] == 'staphylococcus_saprophyticus_coagulase_negative',
                'pathogen'
            ] = 'coagulase_negative_staphylococcus'

        GBD_vs_amr_pathogens = pd.read_csv("FILEPATH")
        df = df.merge(GBD_vs_amr_pathogens, on='pathogen', how='left')
        unlisted_modelling_pathogens = df.loc[df['modelling_side'].isna(), 'pathogen'].unique().tolist()
        print_log_message('These extra pathogens categories are not standard ones found in the above modelling map:')
        print_log_message(unlisted_modelling_pathogens)
        df = df.loc[df['modelling_side'] != 'GBD', ]
        df.drop(columns=['modelling_side'], inplace=True)

        sample_cols = [
            c for c in self.unique_cols if c not in
            [
                'pathogen',
                'pathogen_type',
                'pre_aggregate_pathogen',
                'pathogen_redistributed',
            ]
        ]
        if self.model_type != 'flu_rsv_cw':
            df['can_add_to_poly'] = ~df.pathogen.str.contains(
            	"virus_others|fungi_others|parsite_others|gram_positive_others|gram_negative_others"
            )
            df = df.loc[~(
                df.duplicated(subset=sample_cols, keep=False) & ~df.can_add_to_poly
                & df.groupby(sample_cols).can_add_to_poly.transform(np.any)
            )]
            priority = {
            	'gram_negative_others': 1,
            	'gram_positive_others': 2,
            	'virus_others': 3,
            	'fungi_others': 4,
            	'parasite_others': 5
            }
            rows_to_prioritize = (df.duplicated(subset=sample_cols, keep=False)
                                  & ~df.groupby(sample_cols).can_add_to_poly.transform(np.any))
            sep_df = df.loc[rows_to_prioritize].copy()
            df = df.loc[~rows_to_prioritize]
            sep_df['priority'] = sep_df['pathogen'].map(priority)
            assert sep_df['priority'].notnull().all()
            sep_df = sep_df.sort_values(by='priority').drop_duplicates(
                subset=sample_cols, keep='first'
            )
            sep_df = sep_df.drop('priority', axis='columns')
            df = df.append(sep_df, sort=False)

            df['num_pathogen_per_sample'] = df.groupby(sample_cols, as_index=False)['pre_aggregate_pathogen'].transform('nunique')
            if self.model_type == 'network':
                df['cases'] = df['cases'] / df['num_pathogen_per_sample']
                df['deaths'] = df['deaths'] / df['num_pathogen_per_sample']
                df = df.drop_duplicates(subset=self.unique_cols)    
            elif self.model_type == 'cfr':
                df = df.loc[df['num_pathogen_per_sample'] == 1, ]
            df.drop(columns=['can_add_to_poly', 'pre_aggregate_pathogen', 'num_pathogen_per_sample'], inplace=True)
            self.unique_cols.remove('pre_aggregate_pathogen')
        else:
            df['priority'] = 1
            df.loc[df.pathogen != self.keep_pathogens[0], 'priority'] = 2
            df = df.sort_values(by='priority').drop_duplicates(
                subset=sample_cols, keep='first')

        if self.model_type == 'network':
            others = df.loc[
                (df.source.isin(self.other_node_sources)
                 | df.nid.isin(["NUMBERS"]))
                & (~df.pathogen.isin(self.keep_pathogens)), :].groupby('pathogen')['cases'].agg('sum').reset_index()
            others['syndrome'] = self.infectious_syndrome
            others.to_csv("FILEPATH")
            df.loc[
                (df.source.isin(self.other_node_sources.union({"SOURCE"})) | (
                    df.nid.isin(["NUMBERS"]))
                 ) & (~df.pathogen.isin(self.keep_pathogens)),
                'pathogen'
            ] = 'true_other'
            df = df.loc[df.pathogen.isin(self.keep_pathogens + ['true_other'])]
            df.loc[df.pathogen == 'true_other', 'pathogen'] = 'other'

        elif self.model_type in ['cfr', 'flu_rsv_cw']:
            df.loc[~df.pathogen.isin(self.keep_pathogens), 'pathogen'] = 'other'

        if self.model_type == 'flu_rsv_cw':
            srcs = df.groupby('source').pathogen.apply(
                lambda x: self.keep_pathogens[0] in x.unique()).reset_index()
            df = df.loc[df.source.isin(srcs.loc[srcs.pathogen, 'source'].unique())]

        groupby_cols = [col for col in df.columns if col not in ['deaths', 'cases']]
        df = df.groupby(groupby_cols, as_index=False)['deaths', 'cases'].sum()
        return df

    def prep_age_sex_weights(self, ref_df, pop_df, value_col):
        ref_df = ref_df.loc[ref_df[value_col].notnull()]
        ages = get_cod_ages()
        ref_df = ref_df.loc[
            ref_df.age_group_id.isin(ages.age_group_id.unique().tolist()) &
            ref_df.sex_id.isin([1, 2])
        ]
        if value_col == 'cases':
            age_weight_sources = self.age_weight_sources - self.deaths_only_sources
        else:
            age_weight_sources = self.age_weight_sources
        ref_df = ref_df.loc[ref_df.source.isin(age_weight_sources)]
        missing_ages = set(ages.age_group_id) - set(ref_df.age_group_id)
        if len(missing_ages) > 0:
            warn(
                f"The reference data for age/sex weights is missing the following"
                f" age_group_ids: {missing_ages}"
            )

        group_cols = ['nid', 'location_id', 'year_id',
                      'age_group_id', 'sex_id', 'pathogen']
        ref_df = ref_df.groupby(group_cols, as_index=False)[value_col].sum()
        ref_df = ref_df.append(
            ref_df.groupby(
                [c for c in group_cols if c != 'pathogen'],
                as_index=False
            )[value_col].sum().assign(pathogen='all')
        )
        weight_cols = ['age_group_id', 'sex_id', 'pathogen']
        ref_df_sq = create_square_df(
            ref_df.append(pd.DataFrame(missing_ages, columns=['age_group_id'])),
            weight_cols
        )
        ref_df_sq = ref_df_sq.dropna()
        ref_df_sq = ref_df[['nid', 'location_id', 'year_id']].drop_duplicates()\
            .assign(temp=1)\
            .merge(ref_df_sq.assign(temp=1), on='temp')
        ref_df = ref_df.merge(
            ref_df_sq, how='outer',
            on=['nid', 'location_id', 'year_id', 'age_group_id',
                'sex_id', 'pathogen']
        )
        ref_df[value_col] = ref_df[value_col].fillna(0)

        ref_df = add_population(ref_df, pop_df=pop_df)
        report_if_merge_fail(
            ref_df, 'population',
            ['location_id', 'age_group_id', 'sex_id', 'year_id'])
        ref_df = ref_df.groupby(weight_cols, as_index=False)[
            [value_col, 'population']].sum()
        ref_df['weight'] = ref_df[value_col] / ref_df['population']
        ref_df = ref_df[weight_cols + ['weight']]
        return ref_df

    def add_country_location_id(self, df):
        country_locs = get_country_level_location_id(
            df.location_id.unique().tolist(),
            get_current_location_hierarchy(
                location_set_version_id=self.conf.get_id("location_set_version"),
                **self.cache_kwargs
            )
        )
        df = df.merge(country_locs, how='left', on='location_id', validate='many_to_one')
        report_if_merge_fail(df, 'country_location_id', 'location_id')
        return df

    def age_sex_split(self, df, pop_df, weights, value_col):
        pop_id_cols = ['location_id', 'age_group_id', 'sex_id', 'year_id']

        print_log_message("    Prepping age map")
        age_detail_map = getcache_age_aggregate_to_detail_map(
            **self.cache_kwargs
        )
        sex_detail_map = pd.DataFrame(
            columns=['agg_sex_id', 'sex_id'],
            data=[
                [3, 1],
                [3, 2],
                [9, 1],
                [9, 2],
                [4, 1],
                [4, 2],
                [1, 1],
                [2, 2]
            ]
        )
        detail_maps = {
            'age_group_id': age_detail_map,
            'sex_id': sex_detail_map
        }

        print_log_message("    Prep which distributions should be used")
        pathogen_to_weight_pathogen_map = df[['pathogen']].drop_duplicates()\
            .assign(dist_pathogen=lambda d: d['pathogen'])
        pathogen_to_weight_pathogen_map.loc[
            ~pathogen_to_weight_pathogen_map.pathogen.isin(weights.pathogen),
            'dist_pathogen'
        ] = 'all'
        pathogen_to_weight_pathogen_map['dist_pathogen'].update(
            pathogen_to_weight_pathogen_map['pathogen'].map(
                self.age_weights_use
            )
        )
        self.pathogen_to_weight_pathogen_map = pathogen_to_weight_pathogen_map
        val_to_dist_maps = {
            'pathogen': pathogen_to_weight_pathogen_map
        }
        split_cols = ['age_group_id', 'sex_id']
        split_inform_cols = ['pathogen']
        value_cols = [value_col]
        start_val = df[value_col].sum()
        start_cols = df.columns.tolist()

        df = self.add_country_location_id(df)
        df['orig_location_id'] = df['location_id']
        df['location_id'] = df['country_location_id']

        print_log_message("    Running split")
        df = relative_rate_split(
            df,
            pop_df,
            weights,
            detail_maps,
            split_cols,
            split_inform_cols,
            pop_id_cols,
            value_cols,
            pop_val_name='population',
            val_to_dist_map_dict=val_to_dist_maps,
            verbose=False
        )
        assert np.isclose(start_val, df[value_col].sum())
        return df[start_cols]

    def split_one_metric(self, df, pop_df, weights, value_col):
        keep_cols = self.unique_cols + [value_col]
        df = df[keep_cols].copy()
        df = df.loc[df[value_col].notnull()]
        if value_col == 'cases':
            df = df.loc[~df.source.isin(self.deaths_only_sources)]
        df = self.age_sex_split(
            df, pop_df, weights, value_col
        )
        self.dem_cols = [c for c in self.unique_cols if c != 'sample_id']
        df = df.groupby(self.dem_cols, as_index=False)[value_col].sum()
        return df

    def format_data(self):
        print_log_message(f"Getting data for model type {self.model_type}")
        df = self.get_data()

        print_log_message(f"Selecting specified pathogens")
        df = self.select_pathogens(df)

        print_log_message(f"Generating age/sex splitting weights")
        pop_df = get_pop(
            pop_run_id=self.conf.get_id("pop_run"), **self.cache_kwargs
        )
        self.death_weights = self.prep_age_sex_weights(
            df.copy(), pop_df, value_col='deaths'
        )
        self.case_weights = self.prep_age_sex_weights(
            df.copy(), pop_df, value_col='cases'
        )

        print_log_message("Splitting deaths and cases")
        if self.model_type == 'cfr':
            df = df.loc[df.deaths.notnull() & df.cases.notnull()]
        if self.infectious_syndrome in get_all_related_syndromes(
            'L1_cns_infection',
            self.infsyn,
            level=2
        ):
            gbs_neonates = df.loc[
                (df.pathogen == 'streptococcus_group_b') & (df.age_group_id.isin([2, 3, 42])), 'nid'
            ].unique().tolist()

        if 'blood_stream_infect' in self.infectious_syndrome: 
            df = df.loc[df['year_id'] != 2021]

        deaths_df = self.split_one_metric(
            df, pop_df, self.death_weights, value_col='deaths'
        )
        cases_df = self.split_one_metric(
            df, pop_df, self.case_weights, value_col='cases'
        )

        print_log_message("Recombining data frames")
        df = pd.merge(
            deaths_df, cases_df, how='outer',
            on=self.dem_cols,
            validate='one_to_one'
        )

        print_log_message(
            f"Subsetting to final data needed for {self.model_type}"
        )
        df = df.loc[df.source.isin(self.final_sources)]
        if self.model_type == 'cfr':
            cfrs = df.groupby('source').apply(lambda x: pd.Series(
                {'cfr': sum(x['deaths']) / sum(x['cases'])})).reset_index()
            sources_to_drop = cfrs.loc[cfrs['cfr'].isin([0, 1]), 'source'].tolist()
            df = df.loc[~df['source'].isin(sources_to_drop), :]
            df.loc[df.source.isin(self.icu_only_sources), 'ICU'] = 'ICU_only'
            df.loc[~df.source.isin(self.icu_only_sources), 'ICU'] = 'mixed'
        elif self.model_type in ['network', 'flu_rsv_cw']:
            if self.infectious_syndrome in get_all_related_syndromes(
                'L1_cns_infection',
                self.infsyn,
                level=2
            ):
                df = df.loc[~(
                    ~df.nid.isin(gbs_neonates)
                    & (df.pathogen == 'streptococcus_group_b')
                    & (df.age_group_id.isin([2, 3]))
                )]
                df = df.loc[~((df["source"].isin(["SOURCES"])) 
                    & (df["year_id"].isin(list(range(2003,2009)))))]
            df.loc[df.source.isin(self.icu_only_sources), 'cases'] = np.NaN
            df = df.loc[~(
                df.source.isin(self.icu_only_sources)
                & df.deaths.isnull())]
            assert df.main_diagnosis.notnull().all()
            assert df.infectious_syndrome.notnull().all()
            community = df.infectious_syndrome == df.main_diagnosis
            microbio = df.main_diagnosis == "microbiology"
            df.loc[community & ~microbio, 'hosp'] = 'community'
            df.loc[~community & ~microbio, 'hosp'] = 'hospital'
            if self.infectious_syndrome in get_all_related_syndromes(
                'L1_respiratory_infection',
                self.infsyn,
                level=2
            ):
                flu_rsv_paths = [
                    'influenza_virus', 'respiratory_syncytial_virus', 'poly_flu_rsv'
                ]
                df = df.loc[
                    ~(~community & ~microbio & df.pathogen.isin(flu_rsv_paths))
                ]
                drop = (microbio & df.pathogen.isin(flu_rsv_paths)
                    & (df.hosp == 'hospital'))
                if drop.any():
                    warn(
                        f"Dropping {drop.sum()} rows with pathogen "
                        f"{self.keep_pathogens[0]} and HAI status to CAI"
                    )
                    df = df.loc[~drop]
        return df
