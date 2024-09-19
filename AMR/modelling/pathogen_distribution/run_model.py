"""
Run network meta-analysis to model the distribution of pathogens
causing a given infectious syndrome.
"""
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import sys
import getpass
import yaml
import warnings
import os
from scipy.stats import t
import pickle
from cod_prep.utils import (
    report_duplicates, report_if_merge_fail, create_square_df,
    wait_for_job_ids
)
import itertools
from cod_prep.claude.claude_io import Configurator, makedirs_safely
from mcod_prep.utils.covariates import merge_covariate
from cod_prep.downloaders import (
    add_age_metadata, get_current_location_hierarchy,
    add_location_metadata, get_cod_ages, pretty_print,
    get_pop, add_population, getcache_age_aggregate_to_detail_map,
    get_country_level_location_id, get_ages,
    prep_age_aggregate_to_detail_map
)
from mcod_prep.utils.causes import get_infsyn_hierarchy, get_all_related_syndromes
from cod_prep.utils import print_log_message, warn_if_merge_fail
from amr_prep.utils.pathogen_formatting import PathogenFormatter
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from multiprocessing import Pool
from functools import partial
from pathogen_model_utils import (
    add_country_location_id, add_model_ages, aggregate_data, PathogenCFR
)
sys.path.append(FILEPATH.format(getpass.getuser()))
from netprop.data import Data
from netprop.dorm_model import DormModel
from netprop.model import Model


CONF = Configurator()


class PathogenNetwork(object):
    """Model a pathogen network"""
    blank_pathogens = ['none', 'unknown']
    out_dir = Path(CONF.get_directory("process_data").format(model_step='03a_pathogen'))
    id_cols = ['location_id', 'year_id', 'age_group_id', 'sex_id', 'hosp']

    def __init__(self, model_version, infectious_syndrome, keep_pathogens,
                 covariates, agg_age_group_ids, cfr_use, age_weights_use,
                 year_ids, cfr_ready=None, factor_covs=None, ref_pathogen=None,
                 aggregate_cols=None, unknown=0.5, age_start=None,
                 age_end=None, study_weights=None, gprior_sd=None):
        self.model_version = model_version
        self.infectious_syndrome = infectious_syndrome
        self.agg_age_group_ids = agg_age_group_ids
        self.keep_pathogens = keep_pathogens
        self.ref_pathogen = ref_pathogen
        self.covariates = covariates + ['intercept']
        self.cfr_use = cfr_use or {}
        self.age_weights_use = age_weights_use
        self.year_ids = year_ids
        self.factor_covs = factor_covs or []
        self.cfr_ready = cfr_ready
        self.aggregate_cols = aggregate_cols or []
        self.unknown = unknown
        self.age_start = age_start
        self.age_end = age_end
        self.study_weights = study_weights
        self.gprior_sd = gprior_sd

        self.conf = Configurator()
        self.cache_kwargs = {'force_rerun': True, 'block_rerun': False, 'cache_results': True}
        self.validate_inputs()
        self.model_dir = PathogenNetwork.out_dir / self.infectious_syndrome / self.model_version
        makedirs_safely(str(self.model_dir))

        self.misc_configs = pd.read_csv(
            f"{FILEPATH}/misc_configs.csv"
        ).set_index('entity', verify_integrity=True)['value'].to_dict()

    def validate_inputs(self):
        """Validate inputs"""
        assert self.conf.config_type == 'amr'

        self.infsyn = get_infsyn_hierarchy(infsyn_set_version = 'gram2')
        assert self.infectious_syndrome in self.infsyn.infectious_syndrome.unique()
        self.target_syndromes = get_all_related_syndromes(self.infectious_syndrome, self.infsyn)

        pathogens = pd.read_csv(
            f"{FILEPATH}/pathogen_modelling_map.csv"
        )
        assert set(self.keep_pathogens) <= set(pathogens.pathogen), f"{set(self.keep_pathogens)-set(pathogens.pathogen)} not an acceptable pathogen"
        if self.ref_pathogen is not None:
            assert self.ref_pathogen in pathogens.pathogen.unique()
            assert self.ref_pathogen in self.keep_pathogens

        ages = get_ages(**self.cache_kwargs)
        assert set(self.agg_age_group_ids) <= set(ages.age_group_id)

    def drop_data(self, df):
        # DATA DROPS/OUTLIERS
        if os.path.isfile(f'{FILEPATH}/outliers.csv'):
            print_log_message("Removing outliers...")
            outliers = pd.read_csv(f'{FILEPATH}/outliers.csv')

            outliers = outliers.drop('outlier_note', axis = 1)

            # Merge data and outliers on keys with an indicator
            merged = pd.merge(df, outliers, on=self.id_cols + ['pathogen', 'source'], how='outer', indicator=True, suffixes=('', '_y'))

            # Check for discrepancies in 'cases' column
            diff_cases = merged[merged['_merge'] == 'both']
            if not diff_cases['cases'].equals(diff_cases['cases_y']):
                print_log_message("There are discrepancies in the 'cases' column between data and outliers. CFRs have changed. Outlier re-evaluation required.")

            # Keep only rows that exist only in df
            df = merged.loc[merged['_merge'] == 'left_only',:]

            # Drop the columns that came from 'outliers' and the merge indicator
            df = df.drop(columns=[col for col in df.columns if '_y' in col] + ['_merge'])

        if self.infectious_syndrome == 'L2_blood_stream_infection':
            # pneumocytosis caused by HIV mostly, treated as a contaminant in this iteration of the study
            df = df.loc[df['pathogen'] != 'pneumocystosis_spp',:]
            self.keep_pathogens = [path for path in self.keep_pathogens if path not in ['pneumocystosis_spp']]
        if self.infectious_syndrome == 'L2_urinary_tract_infection':
            # gram positive others to be treated as a contaminant
            df = df.loc[df['pathogen'] != 'gram_positive_others',:]
            self.keep_pathogens = [path for path in self.keep_pathogens if path not in ['gram_positive_others']]
        if self.infectious_syndrome == 'L2_lower_respiratory_infection':
            # candida and adenovirus to be treated as a contaminant
            df = df.loc[~df['pathogen'].isin(['candida_spp', 'adenovirus']),:]
            self.keep_pathogens = [path for path in self.keep_pathogens if path not in ['candida_spp', 'adenovirus']]
        if self.infectious_syndrome == 'L2_meningitis':
            # cryptococcus and toxoplasma caused by HIV mostly, treated as a contaminant in this iteration of the study
            # also drop entamoeba histolytica, high rates of burden from models are not matching what is found in the literature on this
            # pathogen, evidence of bias
            df = df.loc[~df['pathogen'].isin(['cryptococcus_spp', 'toxoplasma_spp', 'entamoeba_histolytica']),:]
            self.keep_pathogens = [path for path in self.keep_pathogens if path not in ['cryptococcus_spp', 'toxoplasma_spp','entamoeba_histolytica']]
            # drop all ICD-coded GBS data for neonates, not aligning with the literature
            df = df.loc[(df['age_group_id'] != 42) | (df['pathogen'] != 'streptococcus_group_b') |
                (~df['data_type'].isin(['mcod', 'linkage', 'hospital'])),:]
        if self.infectious_syndrome == 'L1_bone_joint_infection':
            # 5/3/24 -- all non-pneumo strep in bone joint ICD data is an artificial byproduct of redistribution leading to modelling errors
            df = df.loc[(~df['pathogen'].isin(['streptococcus_group_b', 'streptococcus_group_a', 'streptococcus_others'])) |
                (~df['data_type'].isin(['hospital', 'mcod', 'linkage', 'claims'])),:]
        return df

    def subset_ages(self, df):
        age_group_table = get_ages()
        good_age_group_ids = df.age_group_id.unique().tolist()
        age_detail_map = prep_age_aggregate_to_detail_map(age_group_table, good_age_group_ids)\
            .query("agg_age_group_id in @self.agg_age_group_ids")
        df = df.loc[df['age_group_id'].isin(age_detail_map['age_group_id'].unique()),:]
        return df

    def calculate_se(self, df):
        # Calculate "props" for the purposes of getting uncertainty
        df = df.groupby(
            self.id_cols + ['source', 'nid', 'pathogen', 'data_type'],
            as_index=False
        ).agg(cases=('cases', 'sum'), redist_flag=('redist_flag', 'max')).reset_index()
        df['total_all_paths'] = df.groupby(
            [c for c in self.id_cols if c != 'pathogen'] + ['nid']
        )['cases'].transform(sum)
        # Drop anywhere that was created solely by squaring
        df = df.loc[df.total_all_paths != 0]
        # Drop any 0s as well - these cannot be used in log ratios
        # and earlier sensitivity tests showed that offsetting
        # didn't have a huge effect
        df = df.loc[df.cases != 0]
        df['prop'] = df['cases'] / df['total_all_paths']
        # From binomial distribution
        df['se'] = np.sqrt(df['prop'] * (1 - df['prop']) * df['total_all_paths'])
        assert df.notnull().values.all()
        df = df.drop(['total_all_paths', 'prop'], axis='columns')
        return df

    def set_composites(self, df):
        """
        Set composite pathogens for specific studies.

        """
        # For ICD coded data indicating Gram negative pathogen, represent pathogen entry as a combination of
        # all gram negative bacteria not otherwise covered in that syndrome's ICD codes, in addition to
        # 'Gram negative others'
        if self.infectious_syndrome in ['L2_blood_stream_infection', 'L2_lower_respiratory_infection']:
            bacteria_types = pd.read_csv(f"{FILEPATH}/bacteria_gram_type_map.csv")
            icd_paths = df.loc[df['data_type'].isin(['claims', 'hospital', 'linkage', 'mcod']), 'pathogen'].unique()
            icd_paths = df.loc[df['data_type'].isin(['claims', 'hospital', 'linkage', 'mcod']), 'pathogen'].unique()
            target_gram_negatives = bacteria_types.loc[(bacteria_types['gram_type'] == 'gram_negative') &\
            (~bacteria_types['pathogen'].isin(icd_paths)) & (bacteria_types['pathogen'].isin(self.keep_pathogens)), 'pathogen'].tolist()
            target_gram_negatives = target_gram_negatives + ['gram_negative_others']
            df.loc[df['pathogen'].isin(['gram_negative_unspecified']), 'pathogen'
            ] = '-'.join(target_gram_negatives)
            self.keep_pathogens = [path for path in self.keep_pathogens if path not in ['gram_negative_unspecified']]
        if self.infectious_syndrome == 'skin_infectious':
            df.loc[(df['nid'] == 468273) & (df['pathogen'] == 'other'), 'pathogen']\
                = 'other-acinetobacter_baumanii-enterococcus_faecalis-enterococcus_spp'
        # Set composite pathogen for gbs_meningitis_lit
        elif self.infectious_syndrome == 'L2_meningitis':
            # extraction includes bacterial meningitis only
            df.loc[
                (df.source == 'SOURCE') & (df.pathogen == 'other'),
                'pathogen'
            ] = '-'.join(set(self.keep_pathogens).union({'other'}) - {'streptococcus_group_b', 'non_polio_enteroviruses', 'virus_others'})
        return df

    def get_matches(self, df, ref, case_def, match_cols):
        """
        Match observations for the crosswalk

        """
        assert ref in df[case_def].unique()
        flu_rsv_paths = [
            'flu', 'rsv'
        ]
        assert ref not in flu_rsv_paths
        alts = [x for x in df[case_def].unique() if x != ref]
        assert df[match_cols].notnull().values.all()
        # First match ref to alts
        matches = []
        alts_pair_to_match = alts
        for alt in alts_pair_to_match:
            match = df.query(f"{case_def} == '{ref}'").merge(
                df.query(f"{case_def} == '{alt}'"),
                on=match_cols, validate='one_to_one'
            )
            matches.append(match)
        else:
            combos = list(itertools.combinations(alts, 2))
        for combo in combos:
            match = df.query(f"{case_def} == '{combo[0]}'").merge(
                df.query(f"{case_def} == '{combo[1]}'"),
                on=match_cols, validate='one_to_one'
            )
            matches.append(match)
        df_matched = pd.concat(matches, sort=True)

        # Calculate log ratios and standard errors using delta method
        df_matched['log_ratio'] = np.log(df_matched['cases_y'] / df_matched['cases_x'])
        df_matched['log_ratio_se'] = np.sqrt(
            (df_matched['se_x'] / df_matched['cases_x'])**2
            + (df_matched['se_y'] / df_matched['cases_y'])**2
        )
        return df_matched

    def get_flu_rsv_crosswalks(self, target_cols):
        cw_dir = self.out_dir / "flu_rsv_crosswalk"
        df = pd.read_csv(
            cw_dir / f"{FILEPATH}/data_orig.csv"
        ).append(
            pd.read_csv(cw_dir / f"{FILEPATH}//data_orig.csv")
        )

        # Calculate ratios
        df['log_ratio'] = np.log((1 - df['frac']) / df['frac'])

        # Set pathogens
        flu_paths = ['influenza_virus']
        rsv_paths = ['respiratory_syncytial_virus']
        non_flu_paths = (set(self.keep_pathogens) - set(flu_paths)).union({'other'})
        non_rsv_paths = (set(self.keep_pathogens) - set(rsv_paths)).union({'other'})
        df.loc[df.pathogen == 'flu', 'pathogen_x'] = 'influenza_virus'
        df.loc[df.pathogen == 'rsv', 'pathogen_x'] = 'respiratory_syncytial_virus'
        df.loc[df.pathogen == 'flu', 'pathogen_y'] = '-'.join(non_flu_paths)
        df.loc[df.pathogen == 'rsv', 'pathogen_y'] = '-'.join(non_rsv_paths)

        df['hosp'] = 'community'

        # Calculate log_ratio_se using delta method
        df['log_ratio_se'] = np.sqrt(
            ((df['se'] / df['frac'])**2) * 2
        )

        # Make sure we have the columns we need
        df['cases_x'] = df['frac']
        df['cases_y'] = 1 - df['frac']
        df['se_x'] = df['se']
        df['se_y'] = df['se']
        df['data_type_x'] = df['data_type']
        df['data_type_y'] = df['data_type']
        df = df[target_cols]
        assert df.notnull().values.all()
        return df

    def add_covariates(self, df, covariates=None, add_model_age=True):
        if not covariates:
            covariates = [cov for cov in self.covariates if cov != 'intercept']
        if add_model_age:
            df = add_model_ages(df, self.agg_age_group_ids)
        for cov in covariates:
            if cov not in df:
                if cov == 'hosp_continuous':
                    print_log_message(f"Using unknown = {self.unknown}")
                    df['hosp_continuous'] = df['hosp'].map({
                        'community': 0,
                        'unknown': self.unknown,
                        'hospital': 1
                    })
                elif cov == 'recalc_Hib3_coverage_prop':
                    # slightly reprocessed Hib3 vaccination covariate
                    recalc_hib3 = pd.read_csv(f'{FILEPATH}/reestimated_hib_coverage.csv')
                    df = df.merge(recalc_hib3, on = ['year_id', 'location_id', 'age_group_id'])
                elif cov == 'recalc_PCV3_coverage_prop':
                    # slightly reprocessed PCV3 vaccination covariate
                    recalc_pcv3 = pd.read_csv(f'{FILEPATH}/reestimated_pcv_coverage.csv')
                    df = df.merge(recalc_pcv3, on = ['year_id', 'location_id', 'age_group_id'])
                elif cov == 'cumulative_PCV':
                    # cumulative population level PCV3
                    cumulo_pcv = pd.read_csv(f'{FILEPATH}/cumulative_indirect_pcv_coverage.csv')
                    df = df.merge(cumulo_pcv, on = ['year_id', 'location_id'], validate = 'many_to_one')
                elif cov == 'cumulative_hib':
                    # cumulative population level Hib3
                    cumulo_hib = pd.read_csv(f'{FILEPATH}/cumulative_indirect_hib_coverage.csv')
                    df = df.merge(cumulo_hib, on = ['year_id', 'location_id'], validate = 'many_to_one')
                elif cov == 'inpt_util':
                    # inpatient utilization
                    inpt_util = pd.read_csv(f'{FILEPATH}/inpt_util_AMRest.csv')
                    df = df.merge(inpt_util, on = ['year_id', 'location_id'], validate = 'many_to_one')
                else:
                    df = merge_covariate(
                        df, cov, release_id=self.conf.get_id("release")
                    )
                if cov not in ['recalc_Hib3_coverage_prop', 'recalc_PCV3_coverage_prop']:
                    assert df[cov].notnull().all()
                elif cov == 'recalc_Hib3_coverage_prop':
                    assert df['Hib3_coverage_prop'].notnull().all()
                elif cov == 'recalc_PCV3_coverage_prop':
                    assert df['PCV3_coverage_prop'].notnull().all()
        # convert sex_id into one hot encoded variable
        if 'sex_id' in self.covariates:
            df['sex_id'] = df['sex_id']-1
        return df

    def create_predictions_template(self, pathogens, ages, years):
        """
        Create a predictions template. Should always have all of the id_cols
        """
        lh = get_current_location_hierarchy(
            location_set_version_id=self.conf.get_id('location_set_version'))
        locs = lh.query("level == 3").location_id.unique().tolist()
        if 'sex_id' in self.covariates:
            sexes = [1, 2]
        else:
            sexes = [3]
        if 'hosp' in self.covariates or 'hosp_continuous' in self.covariates:
            hosp = ['community', 'hospital', 'unknown']
        else:
            hosp = ['all']
        index = pd.MultiIndex.from_product(
            [locs, ages, pathogens, years, sexes, hosp],
            names=[
                'location_id', 'age_group_id', 'pathogen',
                'year_id', 'sex_id', 'hosp'
            ]
        )
        square_df = pd.DataFrame(index=index).reset_index()
        square_df = self.add_covariates(square_df, add_model_age=True)
        square_df = square_df.drop('age_group_id', axis='columns')

        # Add CFR
        merge_cols = [c for c in self.cfr_obj.cfr_cols if c != "age_group_id"] + ['agg_age_group_id']
        cfr = self.cfr_obj.cfr.rename(columns={'age_group_id': 'agg_age_group_id'})
        square_df = square_df.merge(
            cfr, how='left', validate='many_to_one', on=merge_cols
        )
        # Fill any thing remaining that failed with the "all pathogens" aggregate
        print_log_message(
            f"Filling {square_df.loc[square_df.cfr.isnull(), merge_cols].drop_duplicates()}")
        square_df = square_df.merge(
            cfr.loc[cfr.pathogen == 'all'].drop('pathogen', axis='columns')[
                [c for c in merge_cols if c != 'pathogen'] + ['cfr']
            ], how='left', validate='many_to_one',
            on=[c for c in merge_cols if c != 'pathogen']
        )
        square_df['cfr'] = square_df['cfr_x'].fillna(square_df['cfr_y'])
        square_df = square_df.drop(['cfr_x', 'cfr_y'], axis='columns')
        if self.cfr_ready:
            report_if_merge_fail(square_df, 'cfr', merge_cols)
        else:
            square_df['cfr'] = square_df['cfr'].fillna(1)

        # Add a bit of metadata for easy plotting
        square_df = add_location_metadata(
            square_df, ['location_name', 'super_region_id'])
        square_df['super_region_name'] = square_df['super_region_id'].map(
            lh.set_index('location_id')['location_name_short'].to_dict()
        )
        return square_df

    def one_hot_encode(self, df):
        # Drop all of the factor covs from the list of covariates
        self.covariates = [cov for cov in self.covariates if cov not in self.factor_covs]
        new_factor_covs = []
        for col in self.factor_covs:
            if df[col].dtype == 'float64':
                df[col] = df[col].astype('int')
            add_cols = pd.get_dummies(
                df[col], prefix=col, prefix_sep='', drop_first=True)
            self.covariates = list(set(self.covariates + add_cols.columns.tolist()))
            new_factor_covs = list(set(new_factor_covs + add_cols.columns.tolist()))
            df = pd.concat([df, add_cols], axis=1)
        self.new_factor_covs = new_factor_covs
        return df

    def set_gaussian_priors(self):
        # Set 0 Gaussian priors on all covariates for every bug
        # axes for priors, 0 - mean coefficient estimate, 1 - prior SD
        if self.gprior_sd is not None:
            self.gpriors = {
                pathogen: {cov: [0, self.gprior_sd] for cov in self.covariates if cov != 'intercept'}
                for pathogen in self.keep_pathogens + ['other']
            }
        else:
            self.gpriors = {
                pathogen: None for pathogen in self.keep_pathogens + ['other']
            }
        
        # emphasize the different distributions by sex for UTI
        if 'sex_id' in self.covariates:
            for pathogen in self.keep_pathogens + ['other']:
                self.gpriors[pathogen]['sex_id'][1] = 0.27
        # minimize the influence of temperature on non-acinetobacter pathogens
        if 'mean_temperature' in self.covariates:
            non_acinetobacter = [bug for bug in self.keep_pathogens + ['other'] if bug != 'acinetobacter_baumannii']
            for pathogen in non_acinetobacter:
                self.gpriors[pathogen]['mean_temperature'][1] = 0.000001
        # minimize the influence of PCV vaccination on non-pneumococcal pathogens
        if 'PCV3_coverage_prop' in self.covariates:
            non_strep = [bug for bug in self.keep_pathogens + ['other'] if bug != 'streptococcus_pneumoniae']
            for pathogen in non_strep:
                self.gpriors[pathogen]['PCV3_coverage_prop'][1] = 0.01
        if 'cumulative_PCV' in self.covariates:
            non_strep = [bug for bug in self.keep_pathogens + ['other'] if bug != 'streptococcus_pneumoniae']
            for pathogen in non_strep:
                self.gpriors[pathogen]['cumulative_PCV'][1] = 0.01
        # minimize the influence of PCV vaccination on non-haemophilus pathogens
        if 'Hib3_coverage_prop' in self.covariates:
            non_hib = [bug for bug in self.keep_pathogens + ['other'] if bug != 'haemophilus_influenzae']
            for pathogen in non_hib:
                self.gpriors[pathogen]['Hib3_coverage_prop'][1] = 0.01
        if 'cumulative_hib' in self.covariates:
            non_hib = [bug for bug in self.keep_pathogens + ['other'] if bug != 'haemophilus_influenzae']
            for pathogen in non_hib:
                self.gpriors[pathogen]['cumulative_hib'][1] = 0.01
        # minimize the influence of PCV vaccination on non-neisseria pathogens
        if 'cv_menafrivac' in self.covariates:
            non_meningo = [bug for bug in self.keep_pathogens + ['other'] if bug != 'neisseria_meningitidis']
            for pathogen in non_meningo:
                self.gpriors[pathogen]['cv_menafrivac'][1] = 0.02
        # minimize the influence of HAQi in bone/joint infections
        if self.infectious_syndrome == 'L1_bone_joint_infection':
            for pathogen in self.keep_pathogens:
                self.gpriors[pathogen]['haqi'][1] = 0.15
        # exceptions as documented in appendix
        elif self.infectious_syndrome == 'L2_lower_respiratory_infection':
            self.gpriors['mycobacterium_others']['haqi'][1] = 0.02
        elif self.infectious_syndrome == 'L2_blood_stream_infection':
            self.gpriors['acinetobacter_baumannii']['haqi'][1] = 0.02
        elif self.infectious_syndrome == 'L2_meningitis':
            self.gpriors['streptococcus_group_b']['haqi'][1] = 0.05
            self.gpriors['non_polio_enteroviruses']['haqi'][1] = 0.05

    def set_uniform_priors(self):
        uprior_covs = ['PCV3_coverage_prop', 'cumulative_PCV', 'Hib3_coverage_prop', 'cumulative_hib', 'cv_menafrivac', 'mean_temperature']
        if any(covariate in self.covariates for covariate in uprior_covs):
            self.upriors = {
               pathogen: {cov: [-np.inf, np.inf] for cov in self.covariates if cov != 'intercept'}
               for pathogen in self.keep_pathogens + ['other']
            }

            # set uniform priors as documented in appendix
            if 'PCV3_coverage_prop' in self.covariates:
                self.upriors['streptococcus_pneumoniae']['PCV3_coverage_prop'] = [-np.inf, 0]
            if 'cumulative_PCV' in self.covariates:
                self.upriors['streptococcus_pneumoniae']['cumulative_PCV'] = [-np.inf, 0]
            if 'Hib3_coverage_prop' in self.covariates:
                self.upriors['haemophilus_influenzae']['Hib3_coverage_prop'] = [-np.inf, 0]
            if 'cumulative_hib' in self.covariates:
                self.upriors['haemophilus_influenzae']['cumulative_hib'] = [-np.inf, 0]
            if 'cv_menafrivac' in self.covariates:
                self.upriors['neisseria_meningitidis']['cv_menafrivac'] = [-np.inf, 0]
            if 'mean_temperature' in self.covariates:
                self.upriors['acinetobacter_baumannii']['mean_temperature'] = [0, np.inf]

    def run_models(self, df_matched, read_model_cache, oosv):
        """Run models with optional leave one country out (LOCO) cross validation"""
        if not read_model_cache:
            uprior_covs = ['PCV3_coverage_prop', 'cumulative_PCV', 'Hib3_coverage_prop', 'cumulative_hib', 'cv_menafrivac', 'mean_temperature']
            if any(covariate in self.covariates for covariate in uprior_covs):
                self.priors = {}
                for pathogen in set(self.gpriors.keys()).union(self.upriors.keys()):
                    self.priors[pathogen] = [self.gpriors.get(pathogen, {}), self.upriors.get(pathogen, {})]

                # Initialize a list to store the rows to save priors
                rows = []
                # Loop over each pathogen in the merged_dict
                for pathogen in self.priors:
                    # Loop over each covariate in the Gprior and Uprior dictionaries
                    for covariate in self.priors[pathogen][0]:
                        # Add a row for the Gprior value
                        rows.append({
                            'pathogen': pathogen,
                            'covariate': covariate,
                            'gprior': self.priors[pathogen][0][covariate],
                            'uprior': self.priors[pathogen][1][covariate]
                        })
                # Convert the list of rows to a DataFrame
                priors_out = pd.DataFrame(rows)
                # Save the DataFrame to a csv
                priors_out.to_csv(self.model_dir / "model_priors.csv", index=False)

                dorm_models = [
                    DormModel(
                        name=pathogen, covs=self.covariates,
                        gprior=prior[0], uprior=prior[1])
                    for pathogen, prior in self.priors.items()
                ]
            else:
                # Initialize a list to store the rows to save priors
                rows = []
                # Loop over each pathogen in the merged_dict
                for pathogen in self.gpriors:
                    # Loop over each covariate in the Gprior and Uprior dictionaries
                    for covariate in self.gpriors[pathogen]:
                        # Add a row for the Gprior value
                        rows.append({
                            'pathogen': pathogen,
                            'covariate': covariate,
                            'gprior': self.gpriors[pathogen][covariate],
                        })
                # Convert the list of rows to a DataFrame
                priors_out = pd.DataFrame(rows)
                # Save the DataFrame to a csv
                priors_out.to_csv(self.model_dir / "model_priors.csv", index=False)

                dorm_models = [
                    DormModel(
                        name=pathogen, covs=self.covariates,
                        gprior=gprior)
                    for pathogen, gprior in self.gpriors.items()
                ]
            with open(self.model_dir / "dorm_models.pkl", 'wb') as file:
                pickle.dump(dorm_models, file)
            worker = f"{self.conf.get_directory('amr_repo')}/"\
                f"{FILEPATH}/model_worker.py"

        if oosv:
            # Reset the index here since we're about to do a lot of
            # column assignment
            df_matched = df_matched.reset_index(drop=True)
            # Add iso3 to distinguish countries
            df_matched = add_location_metadata(df_matched, 'iso3', **self.cache_kwargs)
            assert df_matched.iso3.notnull().all()
            df_matched.to_csv(self.model_dir / "input_data.csv", index=False)
            iso3s = sorted(df_matched.iso3.unique().tolist())
            print(f"Running out-of-sample validation with {len(iso3s)} holdouts")
            data_splits = [
                df_matched.assign(
                   train=lambda d: d['iso3'] != iso3,
                   test=lambda d: d['iso3'] == iso3
                ) for iso3 in iso3s
            ] + [df_matched.assign(train=True, test=True)]
            (self.model_dir / "out_of_sample").mkdir(exist_ok=True)
        else:
            data_splits = df_matched.assign(train=True, test=True, holdout='no_holdout')

        if not read_model_cache:
            print_log_message("Launching workers for modelling...")
            jobs = []

            if oosv:
                for holdout in iso3s + ['no_holdout']:
                    jobname = f"modelworker_{self.model_version}_"\
                        f"{self.infectious_syndrome}_{holdout}"
                    params = [
                        self.model_version, self.infectious_syndrome,
                        holdout, self.ref_pathogen
                    ]
                    jid = submit_mcod(
                        jobname, language='python', worker=worker,
                        cores=10, memory="15G", params=params,
                        runtime="10:00:00", logging=True, queue="long.q",
                        log_base_dir=self.model_dir
                    )
                    jobs.append(jid)
                print_log_message("Waiting...")
                wait_for_job_ids(jobs)
                print_log_message("Jobs complete!")
            else:
                df_matched.to_csv(self.model_dir / "input_data.csv", index=False)
                for holdout in ['no_holdout']:
                    jobname = f"modelworker_{self.model_version}_"\
                        f"{self.infectious_syndrome}_{holdout}"
                    params = [
                        self.model_version, self.infectious_syndrome,
                        holdout, self.ref_pathogen
                    ]
                    jid = submit_mcod(
                        jobname, language='python', worker=worker,
                        cores=10, memory="15G", params=params,
                        runtime="10:00:00", logging=True, queue="long.q",
                        log_base_dir=self.model_dir
                    )
                    jobs.append(jid)
                print_log_message("Waiting...")
                wait_for_job_ids(jobs)
                print_log_message("Jobs complete!")

        print_log_message("Reading cached models...")
        models = {}

        if oosv:
            for holdout in iso3s + ['no_holdout']:
                out_file = self.model_dir / {'no_holdout': ''}.get(
                    holdout, "out_of_sample") / f"model_{holdout}.pkl"
                with open(out_file, 'rb') as file:
                    models[holdout] = pickle.load(file)
            data_splits = dict(zip(iso3s + ['no_holdout'], data_splits))
        else:
            for holdout in ['no_holdout']:
                out_file = self.model_dir / {'no_holdout': ''}.get(
                    holdout, "out_of_sample") / f"model_{holdout}.pkl"
                with open(out_file, 'rb') as file:
                    models[holdout] = pickle.load(file)
        return data_splits, models

    def create_beta_df(self, model):
        buglist = []
        cov_name = []
        betas = []
        for bug in model.dorms:
            betas.extend(model.beta[model.dorm_model_index[bug]].tolist())
            cov_name.extend(model.dorm_models[bug].covs)
            buglist.extend([bug] * len(model.dorm_models[bug].covs))

        betas = pd.DataFrame({'pathogen': buglist, 'cov_names': cov_name, 'beta': betas})

        betas['beta_sd'] = np.sqrt(np.diagonal(model.beta_vcov))

        betas['beta_pval'] = 2 * t.pdf(
            -1 * np.abs(betas['beta'] / betas['beta_sd']),
            model.data.shape[0] - model.beta.shape[0]
        )

        betas.loc[betas['beta_pval'] <= 0.05, 'signif'] = True
        betas.loc[betas['beta_pval'] > 0.05, 'signif'] = False
        return betas

    def get_residuals(self, model, newdata, oosv):
        # Load new data into a new data object
        if oosv:
            newdata = newdata.query("test").reset_index().assign(intercept=1)

        newdata_obj = Data.load(
            newdata,
            obs="log_ratio",
            obs_se="log_ratio_se",
            ref_dorm="pathogen_x",
            alt_dorm="pathogen_y",
            dorm_separator="-"
        )
        # Step 1 - construct dorm_model_mats
        # Model matrices (X) for each pathogen
        dorm_model_mats = {
            name: model.dorm_models[name].get_mat(newdata)
            for name in model.dorms
        }
        # X * beta for each pathogen
        dorm_values = model.get_dorm_values(
            beta=model.beta, dorm_model_mats=dorm_model_mats)
        # Weights w for each pathogen
        ref_dorm_weights = model.get_dorm_weights(newdata_obj.ref_dorm)
        alt_dorm_weights = model.get_dorm_weights(newdata_obj.alt_dorm)
        # Put it all together mathematically
        ref_pred = np.log(np.sum(ref_dorm_weights * dorm_values, axis=1))
        alt_pred = np.log(np.sum(alt_dorm_weights * dorm_values, axis=1))
        newdata['resid'] = newdata_obj.obs.values - (alt_pred - ref_pred)
        return newdata

    def get_prop_residuals(self, df, models, oosv):
        # Some studies have no matches and therefore can't
        # provide useful proportion info
        df = df.loc[df.iso3.isin(models.keys())]
        df = df.assign(
            pathogen_x=self.ref_pathogen,
            pathogen_y=lambda d: d['pathogen'],
            # Log ratio is arbitrary here - it will be added onto
            # the prediction to get the residual and subtracted off below
            # to get the prediction back again
            log_ratio=lambda d: np.log(d['cases']),
            log_ratio_se=1 # arbitrary
        )
        prop_resids = pd.concat([
            self.get_residuals(
                models[holdout],
                df.assign(
                    test=lambda d: d['iso3'] == holdout if holdout != 'no_holdout' else True
                ),
                oosv
            ).assign(holdout=holdout)
            for holdout in models.keys()
        ])
        prop_resids['prop'] = prop_resids['cases'] / prop_resids.groupby(
            self.id_cols + ['nid', 'holdout']
        ).cases.transform(sum)
        prop_resids['cases_pred'] = np.exp(
            prop_resids['log_ratio'] - prop_resids['resid']
        )
        prop_resids['prop_pred'] = prop_resids['cases_pred'] / prop_resids.groupby(
            self.id_cols + ['nid', 'holdout']
        ).cases_pred.transform(sum)
        return prop_resids

    def generate_point_predictions(self, model, preds):
        print_log_message('Generating point predictions...')
        preds = self.one_hot_encode(preds)
        preds['intercept'] = 1

        # generate case point predictions
        # reset index for safety before concatenating columns
        preds = preds.reset_index(drop=True)
        # Model.predict returns columns for each pathogen, sorted
        # alphabetically
        pathogens = sorted(self.keep_pathogens + ['other'])
        preds = pd.concat([
            preds, pd.DataFrame(model.predict(preds), columns=pathogens)
        ], axis=1)
        assert preds.notnull().values.all()
        # Now select the correct column for each row
        preds['prop_cases'] = preds.apply(lambda x: x[x['pathogen']], axis=1)
        preds = preds.drop(pathogens, axis='columns')

        # Calculate prop_deaths
        dem_cols = ['location_id', 'agg_age_group_id', 'sex_id', 'year_id', 'hosp']
        preds['prop_deaths'] = preds['prop_cases'] * preds['cfr']
        preds['prop_deaths'] = preds['prop_deaths'] / preds.groupby(
            dem_cols)['prop_deaths'].transform(sum)

        return preds

    def make_plots(self):
        worker = FILEPATH
        params = [
            self.model_dir, self.infectious_syndrome,
            '--cov_cols'
        ] + self.covariates
        if len(self.factor_covs) > 0:
            params += ['--encode_cols'] + self.new_factor_covs

        # Launch a job
        print_log_message("Launching job for making plots...")
        jobname = f"make_plots_{self.model_version}_{self.infectious_syndrome}"
        jid = submit_mcod(
            jobname, language='r', worker=worker,
            cores=1, memory="10G", params=params,
            runtime="00:30:00", logging=True,
            log_base_dir=self.model_dir
        )
        print_log_message("Waiting...")
        wait_for_job_ids([jid])
        print_log_message("Job complete")

    def run(self, read_data_cache=True, read_model_cache=False, oosv=False):
        #comment out when have cache
        print_log_message("Getting data...")
        if not read_data_cache:
            formatter = PathogenFormatter(
                infsyn_version = 'gram2',
                model_type='network',
                infectious_syndrome=self.infectious_syndrome,
                keep_pathogens=self.keep_pathogens,
                age_weights_use=self.age_weights_use,
                cache_kwargs=self.cache_kwargs
            )
            df = formatter.format_data()

            print_log_message("Saving age/sex weights for vetting")
            pretty_print(formatter.death_weights).to_csv(FILEPATH)
            pretty_print(formatter.case_weights).to_csv(FILEPATH,
                index=False
            )
            df.to_csv(self.model_dir / "raw_data_cache.csv", index=False)
        else:
            print_log_message("Reading from cache")
            df = pd.read_csv(self.model_dir / "raw_data_cache.csv")

        df['nid'] = df['nid'].apply(lambda x: int(x) if type(x) in [int, float] else x)

        print_log_message("Subsetting ages...")
        df = self.subset_ages(df)

        print_log_message("Getting CFRs...")
        self.cfr_obj = PathogenCFR(
            self.infectious_syndrome, self.agg_age_group_ids, self.cfr_use,
            cfr_ready=self.cfr_ready)
        self.cfr_obj.get_cfrs()

        print_log_message("Applying CFRs...")
        df = self.cfr_obj.apply_cfrs(df)
        # Only at this point where we now do not have nulls in the cases columns
        # can we group
        df = df.groupby(
            self.id_cols + ['source', 'nid', 'pathogen', 'data_type'],
            as_index=False).agg(cases=('cases', 'sum'),
            redist_flag=('pathogen_redistributed', 'max')).reset_index()
        print_log_message("Aggregating data...")
        df = aggregate_data(
            df, cols=self.aggregate_cols,
            agg_age_group_ids=self.agg_age_group_ids,
            id_cols=self.id_cols + ['source', 'nid', 'pathogen', 'data_type'])
        print_log_message("Calculating standard errors...")
        df = self.calculate_se(df)

        print_log_message("Applying data drops...")
        df = self.drop_data(df)

        print_log_message("Getting matches...")
        if self.ref_pathogen is None:
            self.ref_pathogen = df.groupby(['pathogen'])['cases'].sum().idxmax()
            print_log_message(
                f"No reference pathogen specified, setting to {self.ref_pathogen}"
            )
        df = self.set_composites(df)
        df_matched = self.get_matches(
            df, self.ref_pathogen, 'pathogen',
            self.id_cols + ['nid', 'source']
        )

        if self.infectious_syndrome == 'L2_lower_respiratory_infection':
            cw = self.get_flu_rsv_crosswalks(df_matched.columns.tolist())
            cw = cw.loc[cw['age_group_id'].isin(df_matched['age_group_id'].unique()), :]
            df_matched = df_matched.append(cw)

        print_log_message("Adding covariates...")
        df = self.add_covariates(df)
        df_matched = self.add_covariates(df_matched)

        if self.infectious_syndrome == 'L2_lower_respiratory_infection':
            # Add Strep pneumo PAFs from vaccine efficacy analysis
            vedata = pd.read_csv(FILEPATH)
            vedata['pathogen_y'] = '-'.join(
                (set(self.keep_pathogens) - {'streptococcus_pneumoniae'}).union({'other'})
            )
            vedata['source'] = vedata['study_type']
            vedata['log_ratio_se'] = vedata['mod_se']

            vedata = vedata.loc[vedata['agg_age_group_id'].isin(df_matched['age_group_id'].unique()), :]

            if 'recalc_PCV3_coverage_prop' in self.covariates:
                vedata = vedata.drop('PCV3_coverage_prop', axis = 'columns')
                recalc_pcv3 = pd.read_csv(FILEPATH)
                vedata = vedata.merge(recalc_pcv3, left_on = ['year_id', 'location_id', 'agg_age_group_id'],
                    right_on = ['year_id', 'location_id', 'age_group_id'])

            if 'recalc_Hib3_coverage_prop' in self.covariates:
                vedata = vedata.drop('Hib3_coverage_prop', axis = 'columns')
                recalc_hib3 = pd.read_csv(FILEPATH)
                vedata = vedata.merge(recalc_hib3, left_on = ['year_id', 'location_id', 'agg_age_group_id'],
                	right_on = ['year_id', 'location_id', 'age_group_id'])

            if 'cumulative_PCV' in self.covariates:
                cumulo_pcv = pd.read_csv(FILEPATH)
                vedata = vedata.merge(cumulo_pcv, on = ['year_id', 'location_id'], validate = 'many_to_one')

            if 'inpt_util' in self.covariates:
                inpt_util = pd.read_csv(FILEPATH)
                vedata = vedata.merge(inpt_util, on = ['year_id', 'location_id'], validate = 'many_to_one')

            if 'cumulative_hib' in self.covariates:
                cumulo_hib = pd.read_csv(FILEPATH)
                vedata = vedata.merge(cumulo_hib, on = ['year_id', 'location_id'], validate = 'many_to_one')


            df_matched = df_matched.append(vedata, sort=False)

            
        print_log_message("Applying study weights...")
        if self.study_weights is not None:
            df_matched['log_ratio_se'] /= df_matched['source'].map(
                self.study_weights).fillna(1)

        print_log_message("Saving input data...")
        df = pretty_print(df, exclude=['nid'])
        if 'sex_id' in self.covariates:
            df.loc[df['sex_id'] == 0, 'sex_label'] = 'Male'
            df.loc[df['sex_id'] == 1, 'sex_label'] = 'Female'
        df['iso3'] = df['ihme_loc_id'].str[0:3]
        df.to_csv(self.model_dir / "data_orig.csv", index=False)

        print_log_message("Creating predictions template...")
        print(self.keep_pathogens)
        preds = self.create_predictions_template(
            self.keep_pathogens + ['other'],
            df.agg_age_group_id.unique().tolist(),
            self.year_ids
        )

        # reframe any recalcualted variables as their base form for modelling
        self.covariates = [c.replace('recalc_', '') for c in self.covariates]

        print_log_message("One-hot encoding factor variables...")
        df_matched = self.one_hot_encode(df_matched)
        df = self.one_hot_encode(df)

        print_log_message("Setting Gaussian priors...")
        self.set_gaussian_priors()

        print_log_message("Setting uniform priors...")
        self.set_uniform_priors()

        if oosv:
            print_log_message("Running model in-sample and with LOCO CV")
        else:
            print_log_message("Running model in-sample")
        data_splits, models = self.run_models(df_matched, read_model_cache, oosv)

        print_log_message("Getting betas and residuals...")
        betas = pd.concat([
            self.create_beta_df(model).assign(holdout=holdout)
            for holdout, model in models.items()
        ])
        # Ratio-space residuals
        if oosv:
            resids = pd.concat([
            self.get_residuals(models[holdout], data_splits[holdout], oosv).assign(
                holdout=holdout
            ) for holdout in models.keys()
            ])
        else:
            resids = self.get_residuals(models['no_holdout'], data_splits, oosv)
        # Proportion-space residuals
        prop_resids = self.get_prop_residuals(df, models, oosv)

        print_log_message("Making predictions...")
        preds = self.generate_point_predictions(models['no_holdout'], preds)

        print_log_message("Saving results")
        resids.to_csv(self.model_dir / "resids.csv", index=False)
        prop_resids.to_csv(self.model_dir / "prop_resids.csv", index=False)
        betas.to_csv(self.model_dir / "betas.csv", index=False)
        preds.to_csv(self.model_dir / "predictions.csv", index=False)

        self.make_plots()
        print_log_message("Done!")


def parse_config(model_version, infectious_syndrome):
    # Read config and parse args
    config_file = pd.read_excel(
        FILEPATH,
        sheet_name='run'
    )
    config = config_file.query(
        f"model_version == '{model_version}' & "
        f"infectious_syndrome == '{infectious_syndrome}'"
    )
    assert len(config) == 1
    config = config.iloc[0].to_dict()
    # Parse lists and dictionaries
    for param in {
        'covariates', 'agg_age_group_ids', 'keep_pathogens',
        'cfr_use', 'age_weights_use', 'year_ids', 'factor_covs',
        'aggregate_cols', 'study_weights'
    }.intersection(set(config.keys())):
        if not config[param] == 'None':
            config[param] = str(config[param]).split(',')
            if param in ['agg_age_group_ids', 'year_ids']:
                if ":" in config[param][0]: 
                    years = str(config[param][0]).split(":")
                    year_start = int(years[0])
                    year_end = int(years[1])
                    years = list(range(year_start,year_end+1))
                    config[param] = [int(x) for x in years]
                else: 
                    config[param] = [int(x) for x in config[param]]
            if param in ['cfr_use', 'age_weights_use', 'study_weights']:
                param_to_value_type = {
                    'cfr_use': str, 'age_weights_use': str,
                    'study_weights': float
                }
                value_type = param_to_value_type[param]
                config[param] = {
                    x.split(':')[0]: value_type(x.split(':')[1]) for x in config[param]
                }
        else:
            config[param] = None
    for param in ['ref_pathogen', 'age_start', 'age_end', 'unknown', 'gprior_sd']:
        if config[param] == 'None':
            config[param] = None

    for param in ['cfr_ready']:
        config[param] = bool(config[param])

    # Deduce keep_pathogens
    if 'keep_pathogens' not in config.keys() or config['keep_pathogens'] is None:
        syn_path = pd.read_csv(FILEPATH)
        syn_path = syn_path.query(
            f"infectious_syndrome == '{config['infectious_syndrome']}'"
        )
        config['keep_pathogens'] = syn_path['pathogens'].iloc[0].split(', ')
    return config


def save_config(out_dir, config):
    with open(f'{out_dir}/config.yml', 'w') as outfile:
        yaml.dump(config, outfile)


if __name__ == '__main__':
    model_version = str(sys.argv[1])
    infectious_syndrome = str(sys.argv[2])
    read_data_cache = str(sys.argv[3])
    read_model_cache = str(sys.argv[4])
    oosv = str(sys.argv[5])
    read_data_cache = read_data_cache == 'True'
    read_model_cache = read_model_cache == 'True'
    oosv = oosv == 'True'
    config = parse_config(model_version, infectious_syndrome)
    print_log_message(
        f"You submitted the following config: {config}"
    )
    network = PathogenNetwork(**config)
    save_config(network.model_dir, config)
    network.run(read_data_cache=read_data_cache, read_model_cache=read_model_cache, oosv=oosv)
