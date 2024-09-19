from http.client import TEMPORARY_REDIRECT
from typing import Any
import os, re, math, datetime
import pandas as pd
import numpy as np
import yaml
import math
import getpass
import json
import pickle
import logging
user_name = getpass.getuser()
from pathlib import Path
from itertools import product
from scipy.special import expit
import plotly.graph_objects as go
import plotly.express as px
import plotly.subplots as sp
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from xspline import XSpline

from cod_prep.claude.claude_io import Configurator
from cod_prep.downloaders.ages import *
from cod_prep.downloaders import get_current_location_hierarchy, add_location_metadata
from mcod_prep.utils.covariates import merge_covariate
from cod_prep.utils import report_duplicates, report_if_merge_fail, create_square_df
from db_tools import ezfuncs

from regmod.variable import Variable
from regmod.models import BinomialModel
from regmod.variable import Variable, SplineVariable
from regmod.utils import SplineSpecs
from regmod.prior import GaussianPrior, UniformPrior, SplineUniformPrior, SplineGaussianPrior
from sklearn.metrics import mean_squared_error, r2_score

class regmod_cfr_generator:
    def __init__(self, syndrome, model_name, limit_result=True, apply_ourtliers=None):
        self.conf = Configurator()
        self.team_dir = Path("ADDRESS")
        self.apply_outliers = apply_ourtliers
        self.replace_in_intercept = pd.DataFrame()
        self.rows_removed = []
        self.fam_outliers = []
        self.syndrome = syndrome
        if 'claims_' in self.syndrome:
            self.syndrome_no_claims = self.syndrome.replace('claims_', '')
        else:
            self.syndrome_no_claims = self.syndrome
        self.model_name = model_name
        self.limited = limit_result
        self.version = datetime.datetime.now().strftime("%m%d")
        if self.apply_outliers:
            self.version = self.version + f"_{self.apply_outliers}outliers"
        self.user_dir = Path(f"ADDRESS")
        self.pathogenpath = "ADDRESS"
        yamlfileloc = 'FILEPATH'
        
        with open(yamlfileloc, 'r') as yaml_file:
            self.params = yaml.safe_load(yaml_file)
        
        self.shared_dir = Path(f"ADDRESS")
        if not os.path.exists(self.shared_dir):
            os.makedirs(self.shared_dir, exist_ok=True)
        self.syndrome_dir = "ADDRESS"
        if not os.path.exists(self.syndrome_dir):
            os.makedirs(self.syndrome_dir, exist_ok=True)
        if not os.path.exists('ADDRESS'):
            os.makedirs('ADDRESS', exist_ok=True)
        if not os.path.exists('ADDRESS'):
            os.makedirs('ADDRESS', exist_ok=True)
        #### SAVE ALL the regmod output models here, all -*- individual, intercept and other/all
        self.regmod_models_dir = 'ADDRESS'
        if not os.path.exists(self.regmod_models_dir):
            os.makedirs(self.regmod_models_dir, exist_ok=True)
        
        self.timestamp = datetime.datetime.now().strftime("%m%d_%H%M")
        log_file_path = "FILEPATH"
        self.log = logging.getLogger("GeneratorLogger")
        self.log.setLevel(logging.DEBUG)
        self.log.propagate = True
        file_handler = logging.FileHandler(log_file_path)
        
        for existing_handler in self.log.handlers[:]:
            self.log.removeHandler(existing_handler)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        self.log.addHandler(file_handler)
        self.log.info(f"================================\t{self.syndrome}  -  {self.model_name}\t================================")
        ### RUN OUTLIER PROCESS IF NEEDED
        if self.apply_outliers:
            if len(self.syndrome)>22:
                sname = self.syndrome[:20]
            else:
                sname = self.syndrome
            if self.model_name != 'all':
                sname = sname+"_"+self.model_name
            try:
                direction = pd.read_excel("FILEPATH", sheet_name=sname)
                if 'direction' not in direction.columns:
                    self.log.warning(f"\t* No outlier comments for {sname} in FILEPATH\n\t\tSKIPPING OUTLIERING")
                    self.direction = pd.DataFrame()
                else:
                    self.direction = direction.dropna(subset=['direction'])
                self.model_pat_dictionary = {}
                for model in direction.model.unique():
                    if model not in ['intercept','individual','family']:
                        continue
                    self.model_pat_dictionary[model] = direction.loc[direction.model==model, 'pathogen'].tolist()
            except:
                self.log.critical(f"**** No previous model outliers directions for {sname} in {self.apply_outliers}_PAT_RANKS\n\t\tSKIPPING OUTLIERING")
                self.direction = pd.DataFrame()
                self.model_pat_dictionary = {}
        if self.limited:
            self.results_dir = Path(f"ADDRESS")
        else:
            self.results_dir = Path(f"ADDRESS")
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir, exist_ok=True)
        self.outlier_dir = Path(f"ADDRESS")
        if not os.path.exists(self.outlier_dir):
            os.makedirs(self.outlier_dir, exist_ok=True)
        
        gbdfilter = pd.read_csv('FILEPATH')
        kp_non_gbd = list(gbdfilter.loc[gbdfilter['modelling_side'].str.contains('AMR'), 'pathogen'])
        droplist = ['staphylococcus_epidermidis_coagulase_negative',
                    'staphylococcus_haemolyticus_coagulase_negative',
                    'staphylococcus_hominis_coagulase_negative',
                    'staphylococcus_others_coagulase_negative']
        if self.syndrome_no_claims == 'meningitis':
            self.keep_pathogens_non_GBD = kp_non_gbd
        elif self.syndrome_no_claims == 'urinary_tract_infection':
            self.keep_pathogens_non_GBD = [p for p in kp_non_gbd if p not in droplist]
        else:
            self.keep_pathogens_non_GBD = [p for p in kp_non_gbd if p not in droplist+['staphylococcus_saprophyticus_coagulase_negative']]
        self.thresholds = {'haqi_cutoff':0.7, 'min_low_haqi_rows':300, 'min_tot_rows':500}
        self.base_params = self.get_dict_from_yaml("all_models")
        self.model_params = self.get_dict_from_yaml(model_name)
        # Get the age conversion dictionary set up
        age_conversion = getcache_age_aggregate_to_detail_map()
        age_conversion = add_age_metadata(age_conversion, ['age_group_name'])
        self.age_conversion_dict = {}
        self.age_names_dict = {42: 'Neonatal',
                               179: 'Post neonatal to 5',
                               1: 'Under 5',
                               239: '5 to 49',
                               25: '50 to 69',
                               26: '70+',
                               }
        self.model_ages = self.get_value('age_categories')
        self.age_color_dict = {age_name: color for age_name, color in zip(list(self.age_names_dict.values()), px.colors.qualitative.D3)}
        for age_id in self.model_ages:
            age_id = int(age_id)
            self.age_conversion_dict[age_id] = age_conversion.loc[age_conversion['agg_age_group_id']==age_id, 'age_group_id'].tolist()
        
        self.sourcelist = self.get_value('source_list')
        self.covariates = self.base_params['always_covariates']
        # Get model types and related covariates
        self.cov_params = self.get_dict_from_yaml("covariates")
        self.covs = ['haqi']
        
        #### GET DATA AND VARIABLES
        self.data = self.get_data()
        self.data['pathogen_type'] = self.data['pathogen_type'].apply(lambda x: re.sub(r'[\\/:"*?<>|]', '_', str(x)))
        # Getting pathogens to estimate
        pathogen_file = Path('FILEPATH')
        estimate_pathogen_dict = pd.read_csv(pathogen_file).set_index('synd')['paths'].str.split(',', expand=False).to_dict()
        self.estimate_pathogen_dict = {key[key.find("_") + 1:]: value for key, value in estimate_pathogen_dict.items()}
        if self.syndrome_no_claims not in self.estimate_pathogen_dict:
            all_unique_pats = set()
            for key, value in self.estimate_pathogen_dict.items():
                all_unique_pats.update(value)
            self.estimate_pathogen_dict[self.syndrome_no_claims] = list(all_unique_pats)
        # Create dictionaries
        self.reference_dict = {}
        self.datadict = {}
        self.model_data = {}
        self.binom_df = {}
        self.vacc_preds = {}
        self.indv_pat_list = []
        # Save_syndrome dictionary is for each model and the respective coefficients
        self.save_syndrome = {}
        self.regmod_binom_model = {}
        self.output = {}
        self.results = {}

        # GENERATE MODEL DATA *-other, family and intercept-*
        list_of_models = list(self.params[self.syndrome_no_claims][self.model_name]['model_types'].keys()) + ['individual']
        for i in range(0, len(list_of_models)):
            m = list_of_models[i]
            self.datadict[m] = self.filter_data(mtype=m,)
            if self.datadict[m].empty:
                self.log.error(f"**** {m} had no data within threshold ****")
                self.datadict[m] = None
            else:
                self.reference_dict[m] = self.get_references(m)

    def get_dict_from_yaml(self, name):
        syndrome = self.syndrome_no_claims
        if 'claims_' in syndrome:
            syndrome = syndrome.replace('claims_', '')
        if name in self.params:
            tempdict = self.params.get(name, {})
        elif name in self.params.get(syndrome, {}):
            tempdict = self.params.get(syndrome, {}).get(name, {})
        else:
            tempdict = self.params.get(syndrome, {}).get(self.model_name, {}).get(name, {})
        # Clean the dictionary to have lists instead of strings and booleans instead of strings
        for key, value in tempdict.items():
            if isinstance(value, str) and ',' in value:
                try:
                    tempdict[key] = [int(item.strip()) for item in value.split(',')]
                except ValueError:
                    tempdict[key] = [item.strip() for item in value.split(',')]
            elif isinstance(value, str):
                if bool(value.strip().lower() == 'true'):
                    tempdict[key] == True
                else:
                    tempdict[key] = value.strip()
            elif isinstance(value, str) and value.isnumeric():
                tempdict[key] = int(value)
        return tempdict

    def get_value(self, key):
        if key in list(self.model_params.keys()):
            return self.model_params[key]
        else:
            return None

    def get_data(self):
        if self.syndrome_no_claims in ['bone_joint_infection','peritoneal_and_intra_abdomen_infection','cardiovascular_infection']:
            filename = f'FILEPATH'
        elif self.syndrome_no_claims in ['myocarditis_pericarditis_carditis','respiratory_non_lower','MMO_encephalitis']:
            filename = f'FILEPATH'
        else:
            filename = f'FILEPATH'
        if 'claims_' in self.syndrome:
            filename = filename.replace('FILEPATH', 'FILEPATH')
        df = pd.read_csv('FILEPATH', index_col=0)
        df['agg_age_id'] = None
        for agg_key, age_vals in self.age_conversion_dict.items():
            df.loc[df.age_group_id.isin(age_vals), 'agg_age_id'] = int(agg_key)        
        df['spec_source'] = df.index # This is the source name
        # ENSURE THAT NO UNIFORM DATA COLUMNS IS INCLUDED
        for covariate_key in self.covariates:
            if df[covariate_key].nunique() == 1:
                self.covariates.remove(covariate_key)
                df.drop(columns=[covariate_key], inplace=True)
                self.log.info(f"**** {covariate_key} has only one value, removed from covariates!!! ****")
        staph_coagulase_neg = ['staphylococcus_epidermidis_coagulase_negative',
                    'staphylococcus_haemolyticus_coagulase_negative',
                    'staphylococcus_hominis_coagulase_negative',
                    'staphylococcus_others_coagulase_negative',
                    'staphylococcus_saprophyticus_coagulase_negative',]
        if any(pathogen in df.pathogen.unique() for pathogen in staph_coagulase_neg):
            self.log.info(f"**** COMBINE ALL coagulase_negative_staphylococcus in data ****")
            df['pathogen'] = np.where(df['pathogen'].isin(staph_coagulase_neg), 'coagulase_negative_staphylococcus', df['pathogen'])
        # Ensure deaths are not greater than cases
        df['deaths'] = np.where(df['deaths'] > df['cases'], df['cases'], df['deaths'])
        # Ensure deaths and cases are not 0
        df = df.loc[df['cases'] != 0]
        if self.syndrome == 'myocarditis_pericarditis_carditis':
            self.log.info(f"**** special treatment for myocarditis_pericarditis_carditis ****")
            df = df.loc[df.agg_age_id != 1]
            self.log.info(f"**** Removed under 5 data for myocarditis_pericarditis_carditis ****")
        #### Exclusion of pathogens by syndrome as suggested by Doctors and AMR team
        exclusion_dict = {'blood_stream_infection': ['pneumocystosis_spp','chikungunya_virus',],
                          'lower_respiratory_infection': ['candida_spp','entamoeba_histolytica'],
                          'meningitis': ['cryptococcus_spp','toxoplasma_spp'],
                          'skin_infection': ['entamoeba_histolytica'],
                          'urinary_tract_infection': ['gram_positive_others'],
                          }
        if self.syndrome in exclusion_dict:
            for pathogen in exclusion_dict[self.syndrome]:
                df = df.loc[df.pathogen != pathogen]
                self.log.info(f"* Excluded {pathogen} from data for {self.syndrome} per reccomendation by AMR team.")
        #### HOTFIX for gram_negative_others
        if self.syndrome in ['blood_stream_infection','lower_respiratory_infection']:
            icd_nid = []
            df = df.loc[~((df['pathogen'] == 'gram_negative_others')&(df['nid'].isin(icd_nid)))]
        ##### APPLY OUTLIERS FROM PREVIOUS MODELS
        if self.apply_outliers:
            removethese, model_pat_dict = self.apply_previous_model_outliers(version=self.apply_outliers)
            if removethese is None:
                self.log.critical(f"****************\n\tNO ROWS TO OUTLIER FROM (SUS...): {self.apply_outliers}\n\t****************")
                return df
            original_cols = df.columns
            merged_df = pd.merge(df, removethese, on=['spec_source', 'pathogen', 'agg_age_id', 'year_id'], 
                     how='left', indicator=True)
            rows_removed = merged_df.loc[merged_df['_merge'] == 'both']
            df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])[original_cols]
            if 'individual' in model_pat_dict:
                self.replace_in_intercept = rows_removed.loc[rows_removed.pathogen.isin(model_pat_dict['individual'])]
            #### SAVE THE ROWS REMOVED
            rows_removed = rows_removed.drop(columns=['_merge'])
            self.log.info(f"****************\n\tRows removed from data due to previous model outliers: {rows_removed.shape[0]}\n\t****************")
            rows_removed[original_cols].to_csv("FILEPATH")
            self.rows_removed = rows_removed
        df = df.reset_index()
        return df

    def apply_previous_model_outliers(self, version='0417'):
        finaloutliers = []
        direction = self.direction.copy()
        if len(direction) == 0:
            self.log.info(f"**** No previous model outliers directions for {self.apply_outliers} ****")
            return None, None

        for model in ['intercept', 'individual','family']:
            if model not in self.model_pat_dictionary:
                self.log.info(f"{model} is not in self.model_pat_dictionary, where it is \n\t{self.model_pat_dictionary}")
                continue
            out_directory = Path(f"ADDRESS")
            file = "FILEPATH"
            if os.path.exists(file):
                self.log.info(f"**** Removing outliers from version: {version} --- {model}")
                outliers = pd.read_csv(file)
                direction = direction[['pathogen','direction','age','out_perc']]
                if direction.empty:
                    self.log.info(f"\t* No outliers to remove for {model} ****")
                for pat, dir, age, perc in direction.values:
                    if model == 'individual':
                        if "all_" in pat:
                            continue
                    if model == 'intercept':
                        if pat not in self.model_pat_dictionary[model]:
                            continue
                    if model == 'family':
                        if 'all_' not in pat:
                            continue
                    if age == 'all':
                        age = self.model_ages
                    else:
                        age = [int(a) for a in age.split(',') if a != '']
                    self.log.info(f"\t* Working on PAT:{pat} for {model} with\n\tDIR:{dir}, AGE:{age} direction and P:{perc} percentile.")
                    if dir == 'low':
                        tempdf = outliers.loc[(outliers.pathogen==pat)&(outliers.agg_age_id.astype(int).isin(age))&(
                            outliers.mean_standardized_res<0)] #&(outliers.observed_cfr<0.05)
                    elif dir == 'high':
                        tempdf = outliers.loc[(outliers.pathogen==pat)&(outliers.agg_age_id.astype(int).isin(age))&(
                            outliers.mean_standardized_res>0)] #&(outliers.observed_cfr>0.34)
                    else:
                        self.log.critical(f"**** Invalid direction {dir} for {pat} in {model} ****")
                    n_drop_keeping = int(tempdf.shape[0] / 0.25 * float(perc))
                    if n_drop_keeping < 3:
                        if len(tempdf) < 3:
                            n_drop_keeping = len(tempdf)
                        else:
                            n_drop_keeping = 3
                    total_rows = len(tempdf)
                    self.log.info(f"Total rows for removal in {pat} in {model} in {dir},{age} {perc}={n_drop_keeping} from 25%={total_rows}")
                    tempdf = tempdf.sort_values(by='mean_standardized_res', ascending=False)
                    tempdf = tempdf.head(n_drop_keeping) if dir == 'high' else tempdf.tail(n_drop_keeping)
                    if tempdf.empty:
                        self.log.info(f"\t* tempdf.empty No outliers to remove for {pat} in {model} ****")
                        continue
                    else:
                        if model == 'family':
                            self.fam_outliers.append(tempdf)
                        else:
                            finaloutliers.append(tempdf)
            else:
                self.log.info(f"**** No previous model outliers for {model} ****")
                continue
        self.log.info("Pass-Fail check 1")
        if len(self.fam_outliers) > 0:
            self.fam_out_df = pd.concat(self.fam_outliers)
        try:
            loadstring = "FILEPATH"
        except:
            loadstring = "FILEPATH"
        oldpath = "FILEPATH"
        self.log.info(f"Attempting to load previous outlier rows from\n{loadstring}")
        if os.path.exists(oldpath):    
            old_final_df = pd.read_csv(oldpath)
            self.log.info(f"\tFOUND {len(old_final_df)} ROWS, and Columns are\n\t{old_final_df.columns}")
            finaloutliers.append(old_final_df)
        else:
            self.log.info(f"{oldpath}")
        self.log.info("Pass-Fail check 2")
        if len(finaloutliers) != 0:
            final_df = pd.concat(finaloutliers)
            self.log.info(f"finaloutliers shape: {final_df.shape}")
            final_df.to_csv(
                "FILEPATH", index=False)
            if len(self.model_pat_dictionary) != 0:
                self.log.info(f"\t****\tTHE MODEL_PAT_DICT IS:\n{self.model_pat_dictionary}")
                with open("FILEPATH", 'w') as json_file:
                    json.dump(self.model_pat_dictionary, json_file)
            self.log.info("Pass-Fail check 3")
            return final_df, self.model_pat_dictionary
        else:
            self.log.info("Pass-Fail check 3")
            return None, None

    def produce_summary_table(self, df):
        pathogen_summary = {}
        for pathogen in df.pathogen.unique():
            minordf = df.loc[df.pathogen==pathogen]
            n_rows = minordf.shape[0]
            n_rows_7 = minordf.loc[minordf.haqi<self.thresholds['haqi_cutoff']].shape[0]
            n_cases_7 = minordf.loc[minordf.haqi<self.thresholds['haqi_cutoff']].cases.sum()
            n_cases = minordf.cases.sum()
            n_deaths = minordf.deaths.sum()
            cfr_all = n_deaths/n_cases
            n_loc = minordf.location_id.nunique()
            min_year = int(minordf.year_id.min())
            max_year = int(minordf.year_id.max())
            mdf = minordf.loc[minordf.cases > 5]
            n_rows_5 = mdf.loc[mdf.haqi<self.thresholds['haqi_cutoff']].shape[0]
            n_rows_5_total = mdf.shape[0]
            n_ages = minordf.agg_age_id.nunique()
            if "ICU" in minordf.columns:
                ICU_cases = minordf.loc[minordf.ICU=="ICU_only"].cases.sum()
                mixed_cases = minordf.loc[minordf.ICU=="mixed"].cases.sum()
                cfr_icu = minordf.loc[minordf.ICU=="ICU_only"].deaths.sum()/ICU_cases
                prop_icu = ICU_cases/n_cases
                pathogen_summary[pathogen] = [n_cases, n_rows, n_rows_7, n_cases_7, n_rows_5, n_rows_5_total, n_deaths, cfr_all,
                                              n_loc, min_year, max_year, n_ages,
                                              ICU_cases, mixed_cases, cfr_icu, prop_icu]
                sumdf = pd.DataFrame.from_dict(pathogen_summary, orient='index',
                                   columns=['n_cases','n_rows','n_rows<0.7haqi', 'n_cases<0.7haqi','n_rows>=5cases<0.7haqi', 'n_rows>=5cases','deaths','cfr_all',
                                            'n_loc', 'first_year', 'last_year', 'n_ages',
                                            'ICU only cases','mixed_cases','cfr_icu','prop_icu'])
            else:
                pathogen_summary[pathogen] = [n_cases, n_rows, n_rows_7, n_cases_7, n_rows_5, n_rows_5_total, n_deaths, cfr_all,
                                              n_loc, min_year, max_year, n_ages,]
                sumdf = pd.DataFrame.from_dict(pathogen_summary, orient='index',
                                    columns=['n_cases','n_rows','n_rows<0.7haqi', 'n_cases<0.7haqi','n_rows>=5cases<0.7haqi', 'n_rows>=5cases','deaths','cfr_all',
                                             'n_loc', 'first_year', 'last_year', 'n_ages',])
        sumdf.sort_values(by='n_cases', ascending=False, inplace=True)
        sumdf['model'] = np.where(sumdf.n_ages == max(sumdf.n_ages), "individual", "intercept")
        sumdf['model'] = np.where(sumdf['n_rows']<=200,"other/all", sumdf['model'])
        self.log.info(f"These pathogens had less than 200 rows of data:\n\t{sumdf.loc[sumdf.n_rows<=200].index}")
        self.log.info(f"These pathogens did had on average less than 0.9 cases per row:\n\t{sumdf.loc[sumdf.n_cases < (sumdf.n_rows*0.9)].index}")
        sumdf['model'] = np.where(sumdf.n_cases < (sumdf.n_rows*0.9), "other/all", sumdf['model'])
        self.top_pathogens = sumdf.index.tolist()[:10]
        tosave = sumdf[[c for c in sumdf.columns if c not in ['n_cases<0.7haqi',
                                                           'mixed_cases','cfr_icu','prop_icu',]
            ]].copy().rename(columns={'n_rows':'Data points', 'n_rows<0.7haqi':'Low (<0.7haqi) datapoints','n_cases':'Total cases',
                                'n_loc':'Unique locations','first_year':'First year','last_year':'Last year', 'deaths':'Total deaths',
                                'cfr_all':'Observed CFR','model':'Model used'
                                })
        tosave.to_csv("FILEPATH")
        return sumdf

    def set_variable(self, column_name, prior_type, params):
        if prior_type == "gaussian":
            if params == 'None':
                p = [GaussianPrior(mean=0.0, sd=0.1)]
            else:
                p = [GaussianPrior(mean=params[0], sd=params[1])]
        elif prior_type == "uniform":
            if params == 'None':
                p = [UniformPrior(lb=-np.inf, ub=np.inf)]
            else:
                p = [UniformPrior(lb=params[0], ub=params[1])]
        else:
            raise ValueError("Invalid type, must be 'gaussian' or 'uniform'.")
        var = Variable(column_name, priors=p)
        return var

    def set_family(self, funcdf, other=False):
        if other:
            funcdf['pathogen'] = 'all'
        else:
            self.pat_family_dict = {pat: f"all_{fam}" for pat, fam in funcdf[['pathogen', 'pathogen_type']].drop_duplicates().values}
            funcdf['pathogen'] = funcdf['pathogen_type'].apply(lambda x: f"all_{x}")
        return funcdf

    def filter_data(self, mtype=None):
        df = self.data.copy()
        if self.apply_outliers:
            if mtype == 'intercept':
                df = pd.concat([df, self.replace_in_intercept])
        if any(df.isnull().sum() > 0):
            self.log.warning(f"\t *-*-* Missing values in data:\n{df.columns[df.isnull().sum() > 0]}")
        if self.get_value("hosp"):
            if self.get_value("hosp") == 'community':
                df = df.loc[df.hosp != 'hospital']
            elif self.get_value("hosp") == 'hospital':
                df = df.loc[df.hosp != 'community']
            else:
                raise ValueError("Invalid value for hosp parameter, must be 'community' or 'hospital'.")
            if 'hosp' not in self.covariates:
                self.covariates.append('hosp')
        df = df[[covariate for covariate in self.covariates if ':' not in covariate] + [
            'deaths','cases','pathogen_type', 'location_id', 'spec_source', 'year_id']]
        self.total_rows = df.shape[0]
        if mtype == 'other':
            self.log.info(f"Number of observations including 'other' >>>> {df.shape[0]}")
            self.log.info(f"Number of observations without 'other' >>>> {df.loc[df.pathogen!='other'].shape[0]}")
            df = self.set_family(df, other=True)
        else:
            df = df.loc[df.pathogen!='other']
            if mtype == 'family':
                df = self.set_family(df, other=False)
                if self.syndrome_no_claims == 'urinary_tract_infection':
                    self.log.info(f"\t* Custom removal for urinary_tract_infection FAMILY, DOCS RECS")
                    df = df.loc[(df['pathogen'] != 'all_virus')]
                if len(self.fam_outliers) > 0:
                    original_cols = df.columns
                    self.log.info(f"\t* Applying previous model outliers for FAMILY from {self.apply_outliers}")
                    self.log.info(f"\t* Removing {len(self.fam_out_df)} rows from {len(df)}...")
                    merged_df = pd.merge(df, self.fam_out_df, on=['spec_source', 'pathogen', 'agg_age_id', 'year_id'], 
                     how='left', indicator=True)
                    df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])[original_cols]
                    self.log.info(f"\t* New number of rows is {len(df)}")
                    self.log.info(f"\t*\n\t**** New columns =\n{df.columns}\n\t**** Old columns =\n{original_cols}")
            else:
                sum_table = self.produce_summary_table(df)
                if mtype == 'intercept':
                    keep_pathogens = list(sum_table.loc[(sum_table['model'] == 'intercept') | (sum_table['model'] == 'individual')].index)
                    self.log.info(f"Pathogens not considered for intercept:\n{[p for p in sum_table.index if p not in keep_pathogens]}")
                    self.log.info(f"Pathogens in intercept model:\n{keep_pathogens}")
                    m_flag = False
                else:
                    keep_pathogens = list(sum_table.loc[sum_table['model'] == 'individual'].index) #self.thresholds['min_tot_rows']].index)
                    self.log.info(f"Unique pathogens for '{mtype}' models: {keep_pathogens}")
                    m_flag = True
                df = df.loc[df.pathogen.isin(keep_pathogens)]
                if ('streptococcus_pneumoniae' in keep_pathogens) or ('haemophilus_influenzae' in keep_pathogens):
                    df = self.merge_vacc_data(df, keep_pathogens)
                    if 'streptococcus_pneumoniae' not in self.vacc_preds and 'streptococcus_pneumoniae' in keep_pathogens:
                        self.get_vacc_model(df, 'streptococcus_pneumoniae')
                    if 'haemophilus_influenzae' not in self.vacc_preds and 'haemophilus_influenzae' in keep_pathogens:
                        self.get_vacc_model(df, 'haemophilus_influenzae')
                if m_flag:
                    for pat in df.pathogen.unique():
                        temp = df.loc[df.pathogen==pat].copy()
                        if pat == 'streptococcus_pneumoniae':
                            pass
                        else:
                            try:
                                temp = temp.drop(columns = ['PCV3_coverage_prop'], axis = 1)
                            except:
                                pass
                        if pat == 'haemophilus_influenzae':
                            pass
                        else:
                            try:
                                temp = temp.drop(columns = ['Hib3_coverage_prop'], axis = 1)
                            except:
                                pass
                        if len(temp.agg_age_id.unique()) == len(self.model_ages):
                            if pat in self.keep_pathogens_non_GBD:
                                self.datadict[pat] = temp
                                self.indv_pat_list = self.indv_pat_list + [pat]
        df = df.drop(columns = ['pathogen_type'], axis = 1)
        return df

    def limit_sources(self, df):
        df.loc[df['source'].isin(self.sourcelist), 'spec_source'] = 'other'
        df.loc[df['location_id'] == 44767, 'location_id'] = 95
        return df

    def merge_vacc_data(self, df, patlist):
        '''This function merges the vaccine coverage data from the vaccine coverage files.
        '''
        df.loc[df['location_id'] == 44767, 'location_id'] = 95
        vacpath = 'ADDRESS'
        if 'streptococcus_pneumoniae' in patlist:
            self.log.info("- - Merged PCV3 coverage data")
            cumulo_pcv = pd.read_csv('FILEPATH')
            df = df.merge(cumulo_pcv, on = ['year_id', 'location_id'], validate = 'many_to_one')
            df = df.rename(columns ={'cumulative_PCV':'PCV3_coverage_prop'})
            if 'pathogen' in df.columns:
                df.loc[df['pathogen'] != 'streptococcus_pneumoniae', 'PCV3_coverage_prop'] = 0
                self.covs = self.covs + ['PCV3_coverage_prop']
            else:
                try:
                    df.loc[df['pathogen_streptococcus_pneumoniae'] != 1, 'PCV3_coverage_prop'] = 0
                except:
                    patcols = [col for col in df.columns if 'pathogen_' in col]
                    mask = (df[patcols] != 0).all(axis=1)
                    df.loc[mask, 'PCV3_coverage_prop'] = 0

        if 'haemophilus_influenzae' in patlist:
            self.log.info("- - Merged HIB coverage data")
            cumulo_hib = pd.read_csv('FILEPATH')
            df = df.merge(cumulo_hib, on = ['year_id', 'location_id'], validate = 'many_to_one')
            df = df.rename(columns ={'cumulative_hib':'Hib3_coverage_prop'})
            if 'pathogen' in df.columns:
                df.loc[df['pathogen'] != 'haemophilus_influenzae', 'Hib3_coverage_prop'] = 0
                self.covs = self.covs + ['Hib3_coverage_prop']
            else:
                try:
                    df.loc[df['pathogen_haemophilus_influenzae'] != 1, 'Hib3_coverage_prop'] = 0
                except:
                    patcols = [col for col in df.columns if 'pathogen_' in col]
                    mask = (df[patcols] != 0).all(axis=1)
                    df.loc[mask, 'Hib3_coverage_prop'] = 0
        return df

    def create_catcont_interaction(self, df, categorical, continuous, m):
        # Attempt to find "reference" parameter, if none found then default to setting the reference as most common category
        if m in ['intercept', 'family', 'other']:
            reference = self.reference_dict[m][categorical]
        elif categorical == 'pathogen':
            reference = m.replace('_binom', '')
        else:
            reference = self.reference_dict['individual'][categorical]
        for col in df.columns:
            if col == f'{categorical}_{reference}':
                continue
            elif categorical in col:
                if col == 'pathogen_type':
                    continue
                df[f"{col}_{continuous}"] = df[col] * df[continuous]
            else:
                continue
        return df

    def get_dummies_columns(self, df, model):
        if "_binom" in model:
            model = model.replace("_binom", "")
        cat_cols = []
        noncat_cols = []
        reference_cols = []
        if model in self.reference_dict.keys():
            model_ref = self.reference_dict[model]
        else:
            model_ref = self.reference_dict['individual']
            model_ref['pathogen'] = model
        for col in df.columns:
            if col == 'spec_source':
                noncat_cols.append(col)
                continue
            if col in self.cov_params and self.cov_params[col]['cov_type']=="category":
                if col in model_ref and model_ref[col] is not None:
                    reference_cols.append(f"{col}_{model_ref[col]}")
                if col not in cat_cols:
                    cat_cols.append(col)
            else:
                if col not in noncat_cols:
                    noncat_cols.append(col)
        df_encoded = pd.get_dummies(df[cat_cols], columns=cat_cols)
        # Remove reference category.
        for col in reference_cols:
            if col in df_encoded.columns:
                df_encoded = df_encoded.drop(columns=col, axis=1)
        # For columns that are getting one-hot encoded, drop column names from self.covs and replace with new column names
        self.covs = [cov for cov in self.covs if cov not in cat_cols]
        for cov in df_encoded.columns:
            if cov not in self.covs:
                self.covs.append(cov)
        df_encoded = pd.concat([df_encoded, df[noncat_cols]], axis=1)
        # Create interaction terms if needed
        if any(':' in c for c in self.covariates) and model != 'intercept':
            for cov in self.covariates:
                if ":" in cov:
                    cat, cont = cov.split(":")
                    df_encoded = self.create_catcont_interaction(df_encoded, cat, cont, model)
        return df_encoded

    def get_references(self, model):
        ref_dict = {}
        data = self.datadict[model].copy().drop(columns=['location_id', 'year_id']) #spec_source
        for cov in self.covariates:
            if cov not in self.cov_params or self.cov_params[cov]['cov_type'] != "category":
                continue
            try:
                reference = self.cov_params[cov]['reference']
            except:
                reference = None
            # Get the model specific reference for other categories
            if reference == "most_common" or reference == "mode" or reference == None:
                mode = data[cov].value_counts()
                if len(mode) == 0:
                    self.log.error(f"**** {cov} has no values in the data: {data[cov].unique()}")
                    continue
                else:
                    reference = mode.idxmax()
            ref_dict[cov] = reference
        return ref_dict

    def prep_dataframe(self, model, data):
        self.log.info(f"====\tRunning prep_dataframe...")
        # Sort pathogens by LEAST common to MOST common, allowing to remove until only most common is left.
        drop_pathogen_ordered = data['pathogen'].value_counts().sort_values().index.tolist()
        pathogens_to_drop = []
        try:
            data = data.drop(columns=['pathogen_type'], axis=1)
        except:
            pass
        for pathogen in data.pathogen.unique():
            dropcols = [col for col in ['pathogen','deaths','cases','spec_source','location_id','year_id'] if col in data.columns]
            tempdf = data.loc[data.pathogen == pathogen].drop(columns= dropcols, axis=1).copy()
            try:
                tempdf = tempdf.drop(columns=['PCV3_coverage_prop'], axis=1)
            except:
                pass
            try:
                tempdf = tempdf.drop(columns=['Hib3_coverage_prop'], axis=1)
            except:
                pass
            columns_with_same_values = tempdf.columns[tempdf.nunique() == 1]
            if len(columns_with_same_values) > 0:
                self.log.warning(f"* - * - {pathogen} Columns with all the same values: {columns_with_same_values}")
                pathogens_to_drop = pathogens_to_drop + [pathogen]
            # CHECK FOR ALL AGES WITHIN PATHOGEN
            if model in ['intercept', 'family', 'other']:
                allages = [a for a in self.model_ages if a != self.reference_dict[model]['agg_age_id']]
            else:
                allages = [a for a in self.model_ages if a != self.reference_dict['individual']['agg_age_id']]
            if not all([age in tempdf.agg_age_id.unique() for age in allages]):
                self.log.warning(f"* - * - {pathogen} Missing age data: {[age for age in allages if age not in tempdf['agg_age_id'].unique()]}")
                if pathogen not in pathogens_to_drop:
                    pathogens_to_drop = pathogens_to_drop + [pathogen]
        self.model_data[model] = data
        try:
            dummy_data = self.get_dummies_columns(data, model)
            prior_vars = self.get_dict_from_yaml('prior_vars')
        except Exception as e:
            self.log.error(f"**** Error with get_dummies_columns or get_dict_from_yaml for {model}:\n{e}\n**** FAILED TO PREP DATA ****")
            return None, None, None
        priorvars = []
        nonpriorvars = []
        priorlist = list(prior_vars.keys())
        for col in dummy_data.columns:
            matchmade = False
            for keyword in priorlist:
                if keyword in col:
                    if col not in priorvars:
                        priorvars.append((col, prior_vars[keyword]))
                    matchmade = True
                    break
            if matchmade == False:
                if col not in nonpriorvars:
                    nonpriorvars.append(col)
        # Remove deaths and cases from potential variables lists
        nonpriorvars = [var for var in nonpriorvars if var not in ['deaths','cases', 'spec_source', 'location_id', 'year_id']]
        priorvars = [(var[0], var[1]) for var in priorvars if var[0] not in ['deaths','cases', 'spec_source', 'location_id', 'year_id']]
        variables = [Variable(x) for x in nonpriorvars]
        for prior in priorvars:
            if prior[1][1] == 'None':
                variables.append(self.set_variable(prior[0], prior[1][0], prior[1][1]))
            else:
                variables.append(self.set_variable(prior[0], prior[1][0], [float(prior[1][1]), float(prior[1][2])]))
        if model not in ['intercept', 'family', 'other']:
            if model != 'streptococcus_pneumoniae':
                    variables = [v for v in variables if 'PCV3' not in v.name]
            if model != 'haemophilus_influenzae':
                variables = [v for v in variables if 'Hib3' not in v.name]
        if any(['coverage_prop' in var.name for var in variables]):
            self.log.info(f"**** {model} Adding negative slope uniform prior to vacc var ****")
            vacc_vars = [var for var in variables if 'coverage_prop' in var.name]
            self.log.info(f"Vaccine variables for {model} model: {vacc_vars}")
            for var in variables:
                if 'coverage_prop' in var.name:
                    var = SplineVariable(var.name,
                                        spline_specs = SplineSpecs(knots=[0,1],
                                                                knots_type = 'rel_freq',
                                                                degree=2,
                                                                l_linear=False,
                                                                r_linear=False,
                                                                include_first_basis=False,
                                                                ),
                                                                priors = [SplineUniformPrior(lb=-np.inf, ub=0, order=1,)]
                                        )
        if model not in ['intercept', 'family']:
            drop_pathogen_ordered = None
        variables = self.gen_prior_age(variables)
        dummy_data['deaths'] = np.where(dummy_data['deaths'] > dummy_data['cases'], dummy_data['cases'], dummy_data['deaths'])
        dummy_data['observed'] = dummy_data['deaths']/dummy_data['cases']
        dummy_data['observed'] = dummy_data['observed'].astype(float)
        if any(dummy_data.isna().sum() > 0):
            self.log.warning(f"**** {model} has NA values in the data: {dummy_data.columns[dummy_data.isna().sum() > 0]} ****")
            dummy_data = dummy_data.dropna()
        return dummy_data, variables, drop_pathogen_ordered
    
    def prep_model_data(self, model):
        self.log.info(f"====\tRunning prep_model_data...")
        if model in self.datadict.keys():
            data = self.datadict[model].copy().drop(columns=['location_id'])
        else:
            data = self.datadict['individual'].copy().drop(columns=['location_id'])
            model = model.replace("_binom", "")
            data = data.loc[data.pathogen==model]
        self.log.info(f"Original Data {data.shape}")
        model_param_dict = self.get_dict_from_yaml('model_types')
        if model in model_param_dict:
            model_covariates = model_param_dict[model]
        else:
            model_covariates = None
        # Add model specific covariates
        if model_covariates:
            additional = model_covariates.strip()
            if "," in additional:
                add_list = additional.split(",")
            else:
                add_list = [additional]
            for item in add_list:
                if item not in self.covariates and item != None:
                    self.covariates = self.covariates + add_list
                else:
                    continue
        self.log.info(f"Getting list of pathogens for {model}...")
        dummy_data, variables, drop_pathogen_ordered = self.prep_dataframe(model, data)
        if any(['coverage_prop' in var.name for var in variables]):
            self.log.info(f"**** {model} Adding negative slope uniform prior to vacc var ****")
            vacc_vars = [var for var in variables if 'coverage_prop' in var.name]
            variables = [var for var in variables if 'coverage_prop' not in var.name]
            self.log.info(f"Vaccine variables for {model} model: {vacc_vars}")
            self.log.info(f"Other variables for {model} model: {variables}")
            for var in vacc_vars:
                try:
                    variables = variables + [SplineVariable(var.name,
                                            spline_specs = SplineSpecs( knots=[0,1],
                                                                        knots_type = 'rel_freq',
                                                                        degree=3,
                                                                        l_linear=False,
                                                                        r_linear=False,
                                                                        include_first_basis=False,
                                                                        ),
                                            priors = [SplineUniformPrior(lb=-np.inf, ub=0, order=1,)]
                                            )]
                except Exception as e:
                    self.log.error(f"**** {e} ****")
        self.log.info(f"DONE WITH DATA PREP FOR {model}")
        return dummy_data, variables, drop_pathogen_ordered

    def get_coefs(self, model, binom_model, vars):
        self.log.info(f"Getting and saving coefficients for {binom_model} after SUCCESSFUL FITTING")
        result = {[x.name for x in vars][i]: binom_model.opt_coefs[i] for i in range(len(vars))}

        if model in ['intercept', 'family', 'other']:
            model_ref = self.reference_dict[model]
        else:
            model_ref = self.reference_dict['individual']
            model_ref['pathogen'] = model
        model_vars = {key: value for key, value in sorted(result.items())}
        model_vars['references'] = model_ref
        if '_binom' not in model:
            mname = model+"_binom"
        else:
            mname = model
        self.save_syndrome[mname] = model_vars

    def produce_model(self, model, binom_df, variables, mod_pat_drop, notused):
        try:
            variables = [v for v in variables if v.name != 'haqi']
            percentiles = [round(n,2) for n in np.percentile(binom_df['haqi'], [20,40,60])]
            min_haqi = min(binom_df['haqi'])
            max_haqi = max(binom_df['haqi'])
            knotlist = [0, percentiles[0], percentiles[2], 1]
            kt = 'abs'
            rl = False
            ll = False
            u_lb = -2
            prior_size = int(len(binom_df.haqi)*0.8)
            thismodelrows = binom_df.shape[0]
            percent = thismodelrows/self.total_rows
            self.log.info(f"Comparing under 0.7 rows {self.total_rows} other/all total and this model has {thismodelrows} rows")
            if model == 'other':
                u_up = 0
                mean = 0
                sd = 1
            else:
                u_up = 0
                mean = 0
                sd = 0.05/percent
                if sd < 1:
                    sd = 1
            if self.syndrome_no_claims == 'eye_infection':
                variables = [v for v in variables if 'agg_age_id_1' not in v.name]
                variables = variables + [Variable('agg_age_id_1', priors=[GaussianPrior(mean=0, sd=.1)])]
            spline_haqi = SplineVariable('haqi',
                                        spline_specs = SplineSpecs(knots=knotlist,
                                                                    knots_type = kt,
                                                                    degree=3,
                                                                    l_linear=ll,
                                                                    r_linear=rl,
                                                                    include_first_basis=False,
                                                                    ),
                                        priors = [SplineUniformPrior(lb=u_lb, ub=u_up, order=1, size=prior_size),
                                                    SplineGaussianPrior(mean=mean, sd=sd)]
                                        )
            variables = variables + [spline_haqi]
            self.log.info(f"ATTEMPT WITH SPLINE - with knots = {knotlist}=[min, 20%, max]\t{ll},{rl}")
            self.log.info(f"\t\tUniform Prior: LB={u_lb}, UB={u_up}, ORDER=1, SIZE={prior_size}")
            self.log.info(f"\t\tGaussian Prior: MEAN=0, SD={sd}")
            binom_model = BinomialModel(y = 'observed',
                        df = binom_df,
                        weights = 'cases',
                        param_specs={'p': {'inv_link': 'expit', 'variables': variables, }}
                        )
        except Exception as e:
            self.log.error(f">>> BinomialModel() CALL DID NOT WORK, ISSUE WITH SPLINE \n\t**** {e}\n")
        try:
            binom_model.fit()
            self.get_coefs(model, binom_model, variables)
            success = True
        except Exception as e:
            self.log.error(f"******************* Original attempt FAIL - {model} *******************")
            self.log.info(f"List of pathogens to remove is:\n\t{mod_pat_drop}")
            success = False
            alt_data = binom_df.copy()
            if mod_pat_drop is None:
                pass
            else:
                for pat in mod_pat_drop[:-1]:
                    if success:
                        break
                    self.log.info(f"Attempting with {pat} removed")
                    alt_data = alt_data.drop(columns=[col for col in alt_data.columns if pat in col], axis=1)
                    variables = [v for v in variables if pat not in v.name]
                    if pat == 'streptococcus_pneumoniae':
                        variables = [v for v in variables if 'PCV3' not in v.name]
                        alt_data = alt_data.drop(columns=['PCV3_coverage_prop'], axis=1)
                    if pat == 'haemophilus_influenzae':
                        variables = [v for v in variables if 'Hib3' not in v.name]
                        alt_data = alt_data.drop(columns=['Hib3_coverage_prop'], axis=1)
                    binom_model = BinomialModel(y = 'observed',
                                        df = alt_data,
                                        weights = 'cases',
                                        param_specs={'p': {'inv_link': 'expit', 'variables': variables, }}
                                        )
                    try:
                        binom_model.fit()
                        self.get_coefs(model, binom_model, variables)
                        success = True
                    except Exception as e:
                        self.log.error(f"**** Still fail after removing {pat} with:\n{e}")
                        tempdf = binom_df.loc[binom_df['pathogen_'+str(pat)]==1].copy()
                        tempdf = tempdf.drop(columns= [col for col in tempdf.columns if 'pathogen' in col], axis=1)
                        columns_with_same_values = tempdf.columns[tempdf.nunique() == 1]
                        if len(columns_with_same_values) > 0:
                            self.log.error(f"**{pat}** Columns with all the same values: {columns_with_same_values}")
                        else:
                            columns_with_same_values = binom_df.dropna(axis=1).columns[binom_df.dropna(axis=1).nunique() == 1]
                            if len(columns_with_same_values) > 0:
                                self.log.error(f"**{pat}** Columns with all the same values: {columns_with_same_values} AFTER dropping nans")
                            else:
                                self.log.warning(f"{pat} seems okay") 
                        pass
            if success == False:
                if model != 'other':
                    self.log.info("***********\nAttempting with age prior after attempts with pathogen removal\n\t\t***********")
                    for possible_col in binom_df:
                        if 'pathogen_' in possible_col:
                            if possible_col not in [var.name for var in variables]:
                                variables = variables + [Variable(possible_col)]
                    prior_vars = [v for v in self.priordist_list if 'agg_age' in v.name]
                    variables = [v for v in variables if 'agg_age' not in v.name]
                    variables = variables + prior_vars
                    self.log.info(f"Variables to go in model:\n{variables}")
                    if len(notused) > 0:
                        self.log.warning(f"Variables NOT present in {model} model: {notused}")
                    alt_data = binom_df.copy()
                    if mod_pat_drop is None:
                        binom_model = BinomialModel(y = 'observed',
                                            df = alt_data,
                                            weights = 'cases',
                                            param_specs={'p': {'inv_link': 'expit', 'variables': variables, }}
                                            )
                        try:
                            binom_model.fit()
                            self.get_coefs(model, binom_model, variables)
                            success = True
                            self.log.info(f"SUCCESS WITH AGE PRIOR!")
                        except:
                            self.log.critical(f"******************* ERROR WITH FITTING - {model} Even with age prior! *******************")
                            showcols = [c for c in binom_df.columns if 'pathogen' in c or 'agg_age' in c]
                            self.log.error(f"VALUE COUNTS OF DATA\n{binom_df[showcols].value_counts()}\n")
                    else:
                        for pat in mod_pat_drop[:-1]:
                            if success:
                                break 
                            self.log.info(f"Attempting with {pat} removed")
                            alt_data = alt_data.drop(columns=[col for col in alt_data.columns if pat in col], axis=1)
                            variables = [v for v in variables if pat not in v.name]
                            if pat == 'streptococcus_pneumoniae':
                                variables = [v for v in variables if 'PCV3' not in v.name]
                                alt_data = alt_data.drop(columns=['PCV3_coverage_prop'], axis=1)
                            if pat == 'haemophilus_influenzae':
                                variables = [v for v in variables if 'Hib3' not in v.name]
                                alt_data = alt_data.drop(columns=['Hib3_coverage_prop'], axis=1)
                            binom_model = BinomialModel(y = 'observed',
                                                df = alt_data,
                                                weights = 'cases',
                                                param_specs={'p': {'inv_link': 'expit', 'variables': variables, }}
                                                )
                            try:
                                binom_model.fit()
                                self.get_coefs(model, binom_model, variables)
                                success = True
                                self.log.info(f"SUCCESS WITH AGE PRIOR!!!!!!!!!!!!!!!!!")
                            except:
                                self.log.critical(f"******************* ERROR WITH FITTING - {model} *******************\n\tEven with age prior!")
                                showcols = [c for c in binom_df.columns if 'pathogen' in c or 'agg_age' in c]
                                self.log.error(f"VALUE COUNTS OF DATA\n{binom_df[showcols].value_counts()}\n")
                else:
                    self.binom_df['other'] = binom_df
                    obs_data = binom_df.copy()
                    obs_data.rename(columns={'observed': 'observed_cfr'}, inplace=True)
                    cols_to_collapse = [c for c in obs_data.columns if 'agg_age_id' in c]
                    obs_data['age_name'] = obs_data.apply(lambda row: row[cols_to_collapse].astype(int).idxmax().replace('agg_age_id_', '') if row[cols_to_collapse].any() else '26',
                        axis=1)
                    fig = px.scatter(obs_data, x='haqi', y='observed_cfr',
                        size='cases', color='age_name',
                        title=f'Observed CFR vs HAQI for',
                        labels={'observed_cfr': 'Observed CFR', 'haqi': 'HAQI'},
                        color_discrete_map= self.age_color_dict, opacity = 0.1,
                        marginal_x='histogram',
                        marginal_y='box',
                        )
                    savepath = 'FILEPATH'
                    fig.to_html(savepath)
                    self.log.critical(f"Graph with data scatter saved to >>> {savepath}")
                    self.log.critical(f"******************* ERROR WITH FITTING 'other' model *******************")
        if success == True:
            self.binom_df[model] = binom_df
            self.regmod_binom_model[model] = binom_model
        return binom_model

    def gen_model(self, model, indv_prior_vars = None):
        self.log.info(f"Binomial model started for {model}...")
        binom_df, variables, mod_pat_drop = self.prep_model_data(model)
        self.log.info(f"Data for {model} prepped...")
        if indv_prior_vars is not None:
            variables = indv_prior_vars
        try:
            notused = [var for var in variables if var.name not in binom_df.columns]
            variables = [var for var in variables if var.name in binom_df.columns]
            if len(notused) > 0:
                self.log.warning(f"Variables NOT present in {model} model: {notused}")
        except:
            self.log.error(f"**** {variables} ****")
            self.log.error(f"**** {binom_df.columns} ****")
        variables = variables + [Variable('intercept')]
        if model in ['intercept', 'family', 'other']:
            agecolsneeded = ["agg_age_id_"+str(a) for a in self.model_ages if a != self.reference_dict[model]['agg_age_id']]
        else:
            agecolsneeded = ["agg_age_id_"+str(a) for a in self.model_ages if a != self.reference_dict['individual']['agg_age_id']]
        if not all([col in binom_df.columns for col in agecolsneeded]):
            self.log.critical(f"\n ****** {[col for col in agecolsneeded if col not in binom_df.columns]} ** MISSING **")
            if model not in ['intercept', 'family', 'other']:
                self.log.info(f"Individual model skipped due to lack of age data {model} <<<<<")
                return None
            else:
                self.log.info(f"continuing with {model}...")
        binom_model = self.produce_model(model, binom_df, variables, mod_pat_drop, notused)
        # SAVE MODEL
        try:
            self.log.info(f"\t==== SAVING {self.syndrome_no_claims} - {model} ====")
            if model not in ['intercept', 'family', 'other']:
                name = model+'_individual'
            else:
                name = model
            with open("FILEPATH", "wb") as f:
                pickle.dump(binom_model, f)
        except Exception as e:
            self.log.info(f"Failed to SAVE MODEL! {name}_regmod_object with pickle.dump- {e}")
        return binom_model

    def gen_prior_age(self, variables):
        """Uniform prior for age groups except neonatal where present"""
        var_list = []
        priordist_list = []
        for var in variables:
            var_name = var.name
            if 'agg_age' in var_name:
                if var_name != 'agg_age_id_42':
                    priordist_list = priordist_list+ [self.set_variable(column_name=var_name,
                                                                    prior_type='uniform',
                                                                    params=[-np.inf, 0])]
                else:
                    priordist_list = priordist_list+ [Variable(var_name)]
            else:
                var_list = var_list + [var]
        var_list = var_list + priordist_list
        return var_list

    def gen_one_model(self, pathogen):
        self.log.info(f"Working on - {pathogen}")
        pat_model = self.gen_model(model = pathogen)
        if pat_model is None:
            self.log.error(f"Failed to generate model for {pathogen}")
        else:
            self.individual_models[pathogen] = pat_model
        return pat_model

    def gen_individual_models(self, ):
        self.individual_models = {}
        for pathogen in self.indv_pat_list:
            pat_mod = pathogen+"_binom"
            pat_model = self.gen_one_model(pathogen)
            if pat_mod in self.save_syndrome.keys():
                self.log.info(f"Done with {self.save_syndrome[pat_mod]['references']['pathogen']}")
            else:
                continue
        self.log.info(f"{self.save_syndrome.keys()}")
        self.log.info(f"{self.individual_models.keys()}")
        return self.individual_models

    def human_readable(self, df_encoded, model):
        df_encoded = self.collapse_col(model+'_binom', df_encoded,'agg_age_id')
        if model in ['intercept', 'family']:
            df_encoded = self.collapse_col(model+'_binom', df_encoded,'pathogen')
        else:
            if model == 'other':
                pathogen = 'all'
            else:
                pathogen = model
            df_encoded['pathogen'] = pathogen
        return df_encoded

    def get_real_v_pred(self):
        if self.syndrome_no_claims == 'unspecified_site_infection':
            csv_file = "FILEPATH"
        else:
            excel_file = "FILEPATH"
            writer = pd.ExcelWriter(excel_file, engine='xlsxwriter')
        output = {}
        individual_df = pd.DataFrame()
        for m_name in self.regmod_binom_model.keys():
            self.log.info(f"Working on real vs predicted for {m_name}...")
            result_df = self.binom_df[m_name].copy()
            result_df['intercept'] = 1
            result_df = self.regmod_binom_model[m_name].predict(result_df)
            result_df['residuals'] = result_df['observed']-result_df['p']
            res_sd = np.std(result_df.residuals)
            result_df['standardized_res'] = (result_df['residuals']/res_sd)*np.sqrt(result_df['cases'])
            result_df.rename(columns={'p': 'predict'}, inplace=True)
            if m_name in ['intercept', 'family', 'other']:
                result_df = self.human_readable(result_df, m_name)
                if self.syndrome_no_claims == 'unspecified_site_infection':
                    result_df.to_csv(csv_file.replace('FILEPATH', 'FILEPATH'), index=False)
                else:
                    result_df.to_excel(writer, sheet_name=m_name, index=False)
            else:
                result_df = self.human_readable(result_df, m_name)
                individual_df = pd.concat([individual_df, result_df], ignore_index=True)
            output[m_name] = result_df
        if not individual_df.empty:
            if self.syndrome_no_claims == 'unspecified_site_infection':
                individual_df.to_csv(csv_file.replace('FILEPATH', f'FILEPATH'), index=False)
            else:
                individual_df.to_excel(writer, sheet_name='individual', index=False)
            output['individual'] = individual_df
        if self.syndrome_no_claims == 'unspecified_site_infection':
            self.log.info(f"Real vs Prediction saved to {csv_file}.")
        else:
            writer.save()
            self.log.info(f"Real vs Prediction saved to {excel_file}.")
        return output

    def collapse_col(self, mod, df, col):
        cols_to_collapse = [c for c in df.columns if col in c]
        if len(cols_to_collapse)>1:
            df[col] = df.apply(lambda row: row[cols_to_collapse].astype(int).idxmax().replace(col+'_', '') if row[cols_to_collapse].any() else self.save_syndrome[mod]['references'][col],
                        axis=1)
        elif len(cols_to_collapse)==0:
            df[col] = self.save_syndrome[mod]['references'][col]
        else:
            val = cols_to_collapse[0]
            df[col] = np.where(df[val]==1,
                            val.replace(col+'_',''),
                            self.save_syndrome[mod]['references'][col])
        try:
            df = df.drop(columns=cols_to_collapse)
        except:
            self.log.info(f"Could not remove {cols_to_collapse} for {mod}...")
        return df

    def make_haqi_range_template(self, model=None):
        # Create template dataframe for all haqi values from 0.20 to 0.98 in increments of 0.01
        haqi_col = pd.Series(np.arange(0.2, 0.99, 0.01))
        results = {}
        if model == None:
            model_list = self.regmod_binom_model.keys()
            self.log.info(f"Working on all models:\n\t\t{model_list}")
        else:
            model_list = [model]
        for mod in model_list:
            self.log.info(f"Working on haqi range template for {mod}...")
            pat_col = self.model_data[mod].pathogen.unique()
            age_col = self.model_data[mod].agg_age_id.unique()
            haqi_col = np.arange(0.1, 0.901, 0.001)
            combinations = list(product(age_col, haqi_col, pat_col))
            template = pd.DataFrame(combinations, columns=['agg_age_id', 'haqi', 'pathogen'])
            if 'ICU' in self.model_data[mod].columns:
                template['ICU']='mixed'
                template['ICU_ICU_only'] = 0
            if 'hosp' in self.model_data[mod].columns:
                template['hosp_'+self.get_value("hosp")] = 1
            template['intercept'] = 1
            template['model'] = mod
            self.log.info(f"Template prediction dataframe for {mod} completed.")
            if mod == 'individual':
                for pathogen in [model for model in self.regmod_binom_model.keys() if model not in ['other', 'intercept', 'family']]:
                    self.log.info(f"Working on individual predictions for {pathogen}.")
                    tempdf = template.loc[template.pathogen == pathogen].copy()
                    tempdf = self.get_dummies_columns(df=tempdf, model=mod)
                    if pathogen == 'streptococcus_pneumoniae':
                        tempdf['PCV3_coverage_prop'] = self.vacc_preds[pathogen]['slope']*tempdf['haqi']
                    if pathogen == 'haemophilus_influenzae':
                        tempdf['Hib3_coverage_prop'] = self.vacc_preds[pathogen]['slope']*tempdf['haqi']
                    tempdf = self.regmod_binom_model[pathogen].predict(tempdf)
                    dropcols = [c for c in tempdf.columns if pathogen in c]
                    tempdf = tempdf.drop(columns=dropcols)
                    tempdf.rename(columns={'p': 'predict'}, inplace=True)
                    tempdf = self.collapse_col('other_binom', tempdf,'agg_age_id')
                    tempdf['pathogen'] = pathogen
                    tempdf['agg_age_id'] = tempdf.agg_age_id.astype(int).replace(self.age_names_dict)
                    keepcols = [c for c in tempdf.columns if 'ICU' not in c or 'hosp' not in c]
                    tempdf = tempdf[keepcols]
                    results[pathogen] = tempdf
                    self.results[(pathogen, mod)] = tempdf
                    self.log.info(f"Individual predictions for {pathogen} DONE.")
            else:
                self.log.info(f"\tCalculating predictions for {mod}.")
                template = self.get_dummies_columns(df=template, model=mod)
                if 'streptococcus_pneumoniae' in pat_col:
                    self.log.info(f"Calculating PCV3 coverage for {mod}.")
                    try:
                        template['PCV3_coverage_prop'] = np.where(template['pathogen_streptococcus_pneumoniae']==1,self.vacc_preds['streptococcus_pneumoniae']['slope']*template['haqi'],0)
                    except:
                        allpatcols = [col for col in template.columns if 'pathogen' in col]
                        template['PCV3_coverage_prop'] = np.where((template[allpatcols] == 0).all(axis=1), self.vacc_preds['streptococcus_pneumoniae']['slope'] * template['haqi'], 0)
                if 'haemophilus_influenzae' in pat_col:
                    self.log.info(f"Calculating Hib3 coverage for {mod}.")
                    try:
                        template['Hib3_coverage_prop'] = np.where(template['pathogen_haemophilus_influenzae']==1,self.vacc_preds['haemophilus_influenzae']['slope']*template['haqi'],0)
                    except:
                        allpatcols = [col for col in template.columns if 'pathogen' in col]
                        template['Hib3_coverage_prop'] = np.where((template[allpatcols] == 0).all(axis=1), self.vacc_preds['haemophilus_influenzae']['slope'] * template['haqi'], 0)
                res_df = self.regmod_binom_model[mod].predict(template)
                res_df.rename(columns={'p': 'predict'}, inplace=True)
                dropcols = [c for c in res_df.columns if '_haqi' in c]
                res_df = res_df.drop(columns=dropcols)

                res_df = self.collapse_col(mod+'_binom', res_df,'agg_age_id')
                res_df = self.collapse_col(mod+'_binom', res_df,'pathogen')
                res_df['agg_age_id'] = res_df.agg_age_id.astype(int).replace(self.age_names_dict)
                keepcols = [c for c in res_df.columns if 'ICU' not in c or 'hosp' not in c]
                res_df = res_df[keepcols]
                results[mod] = res_df
                for pathogen in res_df.pathogen.unique():
                    tempdf = res_df.loc[res_df.pathogen == pathogen].copy()
                    self.results[(pathogen, mod)] = tempdf
                self.log.info(f"Predictions for {mod} DONE.")
        return results

    def box_plot_by_pathogen_ordered(self, model=None, group='pathogen', colx='spec_source', coly='standardized_res'):
        """
        Create box plots for each unique value in the 'pathogen' column,
        ordered by the average value of the y-column.
        Parameters:
        - model: model name (default is "other_binom" model)
        - colx: str, the categorical variable column name
        - coly: str, the continuous variable column name
        - save: bool, whether to save the plot or not (default is False)
        - save_path: str, the path where the plot should be saved (default is ".")
        Returns:
        - None (displays box plots/ saves to residual_box_plots folder)
        """
        if model is None:
            data = self.output['other'].copy()
            model = 'other'
        else:
            data = self.output[model].copy()

        if group not in data.columns:
            self.log.warning(f"The {group} column is not present in the DataFrame.")
            if 'ICU_ICU_only' in data.columns:
                group = 'ICU_ICU_only'
            elif 'hosp' in data.columns:
                group = 'hosp'
            else:
                group = 'intercept'
        if 'ICU_ICU_only' in data.columns:
            self.log.info(f"SHOWING {data.ICU_ICU_only.unique()}")
        avg_y_values = data.groupby(group)[coly].mean().sort_values(ascending=False)
        sorted_groups = avg_y_values.index
        fig = px.box(data, x=colx, y=coly, color=group, category_orders={group: sorted_groups})
        fig.update_layout(font=dict(size=20),showlegend=True,
            title_text=f'{self.model_name}: {coly} grouped by {group} and {colx}',
                        shapes=[
                            dict(type='line', y0=0, y1=0, xref='paper', x0=0, x1=1, line=dict(color='black', width=1))
                        ]
        )
        save_folder = "ADDRESS"
        os.makedirs(save_folder, exist_ok=True)
        file_path = "FILEPATH"
        fig.write_html(str(file_path))

    def create_interactive_age_graphs(self, preds):
        for pathogen in preds.pathogen.unique():
            tempdf = preds.loc[preds.pathogen == pathogen].copy()
            unique_models = tempdf.model.unique()
            if len(unique_models) > 1 and 'intercept' in unique_models:
                tempvalue = abs(tempdf.loc[(tempdf.model == 'individual')&(tempdf.agg_age_id == '70+')
                                            ].predict.max()-tempdf.loc[(tempdf.model == 'individual')&(tempdf.agg_age_id == '70+')
                                                                        ].predict.min())
                if  tempvalue < 0.01:
                    excludemod = 'individual'
                else:
                    excludemod = 'intercept'
                if pathogen in self.pathogen_intercept_winner:
                    excludemod = 'individual'
                excluded_index = preds.loc[(preds.pathogen == pathogen) & (preds.model == excludemod)].index
                preds = preds.drop(index = excluded_index)
        unique_pathogens = preds['pathogen'].unique()
        num_pathogens = len(unique_pathogens)
        num_rows = math.isqrt(num_pathogens)
        num_cols = math.ceil(num_pathogens / num_rows)

        fig = sp.make_subplots(rows=num_rows, cols=num_cols, subplot_titles=unique_pathogens)

        for i, pathogen in enumerate(unique_pathogens, start=1):
            pathogen_data = preds[preds['pathogen'] == pathogen]
            row_index = (i - 1) // num_cols + 1
            col_index = (i - 1) % num_cols + 1
            line_fig = px.line(
                pathogen_data,
                x='haqi',
                y='predict',
                color='agg_age_id',
                color_discrete_map=self.age_color_dict,
                labels={'predict': 'Case Fatality Rate', 'haqi': 'HAQI'},
                title=f'{pathogen}'
            )
            for trace in line_fig.data:
                fig.add_trace(trace, row=row_index, col=col_index)
        counter = 0
        for trace in fig.data:
            if counter>=len(preds.agg_age_id.unique()):
                trace.update(showlegend=False)
            counter += 1
        fig.update_layout(font=dict(size=20),#showlegend=True,
                          title_text=f"Case Fatality Rate: HAQI for {self.syndrome}_{self.model_name}",
                          updatemenus=[{
                                    'buttons': [
                                        {'args': ['hovermode', 'x'], 'label': 'Compare Data', 'method': 'relayout'},
                                        {'args': ['hovermode', 'closest'], 'label': 'Hover Closest', 'method': 'relayout'},
                                    ],
                                    'direction': 'down',
                                    'showactive': True,
                                    }],
                        annotations=[dict(font=dict(size=10))]
                        )
        output_folder = "ADDRESS"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        output_file = "FILEPATH"
        fig.write_html(output_file)
        return output_file

    def create_individual_interactive_graphs(self, preds, allmodel=False):
        if allmodel:
            preds = preds.loc[preds['model'].isin(['other'])]
        else:
            preds = preds.loc[preds['model'].isin(['individual', 'intercept'])]
            # If no individual or intercept, do the family model instead
            if preds.empty:
                preds = preds.loc[preds['model'].isin(['family'])]
        self.log.info(f"Working on {preds.model.unique()} models")
        cols = [c for c in self.data.columns if c in ['spec_source', 'pathogen', 'deaths', 'cases', 'haqi', 'agg_age_id', 'ICU']]
        obs_data = self.data[cols]
        obs_data.loc[obs_data['deaths'] > obs_data['cases'], 'deaths'] = obs_data['cases']
        obs_data['observed_cfr'] = obs_data['deaths'] / obs_data['cases']
        
        def min_max_scaling(series, min_bound=4, max_bound=20):
            max_value = series.max()
            min_value = series.min()
            scaled_series = ((series - min_value) / (max_value - min_value)) * (max_bound - min_bound) + min_bound
            return scaled_series
        obs_data['size'] = obs_data.groupby('agg_age_id',group_keys=False)['cases'].apply(lambda x: min_max_scaling(x, 2, 20))
        
        obs_data['age_name'] = obs_data['agg_age_id'].astype(int).replace(self.age_names_dict)
        obs_data = obs_data.drop(columns=['deaths','agg_age_id'])
        if allmodel:
            obs_data['pathogen'] = 'all'
            set_pathogens = ['all']
        else:
            set_pathogens = set(obs_data.pathogen.unique()) & set(preds.pathogen.unique())
        if len(set_pathogens) == 0:
            self.log.critical(f"NO PATHOGENS IN COMMON BETWEEN OBSERVED AND PREDICTED DATA")
        else:
            set_pathogens = [p for p in list(set_pathogens) if p in self.keep_pathogens_non_GBD]
            self.log.info(f"Working on {set_pathogens}")
        
        if allmodel:
            output_path = "ADDRESS"
        else:
            output_path = "ADDRESS"
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        
        hovercols = [c for c in obs_data.columns if c in ['haqi', 'observed_cfr', 'cases', 'spec_source', 'ICU']]
        if 'ICU' in hovercols:
            s_col = 'ICU'
        else:
            s_col = None
        # Loop through each pathogen
        for pathogen in set_pathogens:
            pathogen_preds = preds.loc[preds.pathogen == pathogen]
            obs_data_pathogen = obs_data.loc[obs_data.pathogen == pathogen]
            fig = px.scatter(obs_data_pathogen, x='haqi', y='observed_cfr',
                            size='size', color='age_name', symbol=s_col,
                            title=f'Observed CFR vs HAQI for {self.model_name}_{pathogen}',
                            labels={'observed_cfr': 'Observed CFR', 'haqi': 'HAQI'},
                            hover_data=hovercols,
                            color_discrete_map=self.age_color_dict, opacity = 0.2,
                            marginal_x='histogram',
                            marginal_y='box',
                            )
            # Fix the image paramters
            def update_marker_shape(trace):
                if "ICU_only" in trace.name:
                    trace.update(showlegend=False)
                    trace.update(opacity=0.0)
                    if trace.type == 'scatter':
                        trace.marker.symbol= 'diamond'
                else:
                    if trace.type == 'scatter':
                        trace.marker.symbol= 'circle'
            fig.for_each_trace(update_marker_shape)

            # Add lines for predicted CFR
            for age, age_data in pathogen_preds.groupby('agg_age_id'):
                for model_type, model_data in age_data.groupby('model'):
                    dash_style = 'solid' if model_type == 'intercept' else 'dash'
                    fig.add_scatter(x=model_data['haqi'], y=model_data['predict'],
                                    mode='lines', line=dict(dash=dash_style, color=self.age_color_dict[age]),
                                    name=f'cfr for {age} {model_type}',)
            
            fig.update_layout(font=dict(size=20),updatemenus=[
                    {'buttons': [
                        {'args': ['hovermode', 'x'], 'label': 'Compare Data', 'method': 'relayout'},
                        {'args': ['hovermode', 'closest'], 'label': 'Hover Closest', 'method': 'relayout'},
                    ],
                    'direction': 'down',
                    'showactive': True,}
                ],
                yaxis=dict(range=[-0.1, 1.1]),
            )
            pathogen = re.sub(r'[\\/:"*?<>|]', '_', pathogen)
            graph_filename = "FILEPATH"
            fig.write_html(graph_filename)
        return output_path

    def create_heatmap_graph(self, preds, pathogen=False):
        if pathogen=='all':
            preds = preds.loc[preds['model'].isin(['other'])]
        elif pathogen != False:
            preds = preds.loc[preds['pathogen'] == pathogen]
        else:
            preds = preds.loc[preds['model'].isin(['individual', 'intercept'])]
            if preds.empty:
                return "NO PATHOGENS IN 'INDIVIDUAL' OR 'INTERCEPT' MODELS"
        cols = [c for c in self.data.columns if c in ['spec_source', 'pathogen', 'deaths', 'cases', 'haqi', 'agg_age_id', 'ICU']]
        obs_data = self.data[cols]
        if pathogen == 'all':
            obs_data['pathogen'] = 'all'
        elif pathogen != False:
            obs_data = obs_data.loc[obs_data['pathogen'] == pathogen]
        else:
            obs_data = obs_data.loc[obs_data['pathogen'].isin(preds.pathogen.unique())]
        obs_data.loc[obs_data['deaths'] > obs_data['cases'], 'deaths'] = obs_data['cases']
        obs_data['observed_cfr'] = obs_data['deaths'] / obs_data['cases']
        step = 0.05
        haqi_bins = np.arange(0,1+step, step)
        z, x = np.histogram(
            obs_data['cases'],
            bins=[haqi_bins],
        )
        x = x[:-1] + step / 2
        fig = go.Figure(data=go.Heatmap(z=z,x=x, y=[1],
                                        colorscale='Purples',))

        for pathogen in preds.pathogen.unique():
            pathogen_preds = preds.loc[preds.pathogen == pathogen]
            for age, age_data in pathogen_preds.groupby('agg_age_id'):
                for model_type, model_data in age_data.groupby('model'):
                    dash_style = 'solid' if model_type in ['individual', 'other'] else 'dash'
                    fig.add_scatter(x=model_data['haqi'], y=model_data['predict'],
                                    mode='lines', line=dict(dash=dash_style, color=self.age_color_dict[age]),
                                    name=f'cfr for {age} {model_type}',)
            
            pathogen = re.sub(r'[\\/:"*?<>|]', '_', pathogen)
            fig.write_html("FILEPATH")

    def create_interactive_pathogen_graphs(self, preds, age='5 to 49'):
        # Define linetypes based on the 'model' column
        agecol = 'agg_age_id'
        if agecol not in preds.columns:
            agecol = 'age'

        # Get unique values in agecol
        ipg_df = preds.loc[preds[agecol] == age].copy()
        for pat in ipg_df.pathogen.unique():
            tempdf = ipg_df.loc[ipg_df.pathogen == pat]
            unique_models = tempdf.model.unique()
            if len(unique_models) > 1 and 'intercept' in unique_models:
                excludemod = 'intercept'
                excluded_index = ipg_df.loc[(ipg_df.pathogen == pat) & (ipg_df.model == excludemod)].index
                ipg_df = ipg_df.drop(index = excluded_index)
        n_pats = len(ipg_df.loc[ipg_df['model'].isin(['individual', 'intercept', 'other'])].pathogen.unique())
        
        if n_pats > len(self.top_pathogens):
            self.log.warning(f"Over 15 pathogens to plot - {n_pats}")
            splitlist = [self.top_pathogens, [ipg_df.loc[~ipg_df['pathogen'].isin(self.top_pathogens)]]]
        else:
            splitlist = [ipg_df.pathogen.unique()]
        for n, list in enumerate(splitlist):
            ipg_df = ipg_df.loc[ipg_df['pathogen'].isin(list)]
            fig = px.line(
                ipg_df,
                x='haqi',
                y='predict',
                color='pathogen',
                line_dash='model',
                labels={'predict': 'Case Fatality Rate', 'haqi': 'HAQI'},
            )

            # Update layout
            fig.update_layout(font=dict(size=20),
                xaxis_title='HAQi',
                yaxis_title='Case Fatality Rate',
                title={'text': f'Case Fatality Rate: HAQI for {self.syndrome}_{self.model_name} pathogens',
                    'y': 0.95, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'},
                legend=dict(
                    orientation='v',
                    y=0.5,
                    x=1.02,
                ),
                updatemenus=[
                    {'buttons': [
                        {'args': ['hovermode', 'x'], 'label': 'Compare Data', 'method': 'relayout'},
                        {'args': ['hovermode', 'closest'], 'label': 'Hover Closest', 'method': 'relayout'},
                    ],
                    'direction': 'down',
                    'showactive': True,}
                ],
            )
            if n == 0:
                output_file = "FILEPATH"
                filelist = [output_file]
            else:
                output_file = "FILEPATH"
                filelist.append(output_file)
            fig.write_html(output_file)
        return filelist

    def r_plots_data(self):
        savepath = "FILEPATH"
        finaldf = pd.DataFrame()
        for key, df in self.results.items():
            pat, m = key
            if m not in ['intercept', 'family', 'other']:
                m = 'individual'
            df['model'] = m
            df['pathogen'] = pat
            if 'hosp' in df.columns:
                df = df.loc[df['hosp'] == 'unknown']
            finaldf = pd.concat([finaldf, df], ignore_index=True)
            self.log.info(f"Added {pat} from {m} to finaldf - {finaldf.shape}")
        self.log.info(f"Finaldf finalized, saving to {savepath}")
        finaldf = finaldf.reset_index(drop=True)
        if finaldf.empty:
            self.log.critical(f">>>>>>>>>>>>>>>>>> {self.syndrome} {self.model_name} has no data to plot.")
            raise ValueError(f">>>>>>>>>>>>>>>>>> {self.syndrome} {self.model_name} has no data to plot.")
        finaldf.to_csv(savepath, index=False)
        self.log.info(f"- - R script plot data saved.")
        return finaldf

    def gen_output(self):
        for m in self.params[self.syndrome_no_claims][self.model_name]['model_types'].keys():
            if self.datadict[m] is None:
                continue
            self.log.info(f"\t\t==== Generating {m} model ====")
            self.gen_model(model=m)
            if m == 'other':
                sd = 1
                priordist_list = []
                for var_name in self.save_syndrome['other_binom'].keys():
                    if 'agg_age_id' in var_name:
                        if var_name != 'agg_age_id_42':
                            self.log.warning(f"self.save_syndrome['other_binom'][{var_name}] = {self.save_syndrome['other_binom'][var_name]}")
                            priordist_list = priordist_list+ [Variable(var_name,
                                                                        priors = [GaussianPrior(mean=self.save_syndrome['other_binom'][var_name], sd=sd),
                                                                                    UniformPrior(lb=-np.inf, ub=0)])]
                        else:
                            self.log.warning(f"self.save_syndrome['other_binom'][{var_name}] = {self.save_syndrome['other_binom'][var_name]}")
                            priordist_list = priordist_list+ [Variable(var_name,
                                                                        priors = [GaussianPrior(mean=self.save_syndrome['other_binom'][var_name], sd=sd)])]
                self.priordist_list = priordist_list
            self.log.info(f"\t\t======== {m} model FINISHED ========\n")

        self.log.info(f"\t\t======== Generating individual models ========") 
        self.gen_individual_models()
        self.log.info(f"\t\t======== Individual models FINISHED ========\n")

        try:
            self.log.info(f"\t\t ==== SAVING OBJECT {self.syndrome}_{self.model_name} ====")
            objectspath = "ADDRESS"
            if not os.path.exists(objectspath):
                os.makedirs(objectspath, exist_ok=True)
            with open("FILEPATH", "wb") as f:
                pickle.dump(self, f)
        except Exception as e:
            self.log.info(f"Failed to SAVE OBJECT! With pickle.dump - {e}")

        self.log.info(f"\t\tFunction make_haqi_range_template() run for all models!")
        self.make_haqi_range_template()
        self.log.info(f"Function make_haqi_range_template() FINISHED\n")
        self.log.info(f"\n\t\tFunction get_real_v_pred() run for all models!")
        self.output = self.get_real_v_pred()
        self.log.info(f"Function get_real_v_pred() FINISHED\n")
        
        self.log.info(f"Determining if pathogen should be predicted in individual or intercept models...")
        self.pathogen_intercept_winner = self.pick_model()
        # FIXING AUTOMATED PROCESS FOR SOME PATHOGENS
        for model in self.pathogen_intercept_winner:
            if 'influenza_virus' not in self.pathogen_intercept_winner and self.syndrome_no_claims == 'lower_respiratory_infection':
                self.pathogen_intercept_winner.append('influenza_virus')
            elif 'mycobacterium_others' not in self.pathogen_intercept_winner and self.syndrome_no_claims == 'lower_respiratory_infection' and self.model_name == 'community':
                self.pathogen_intercept_winner.append('mycobacterium_others')
            elif 'enterobacter_spp' not in self.pathogen_intercept_winner and self.syndrome_no_claims == 'lower_respiratory_infection':
                self.pathogen_intercept_winner.append('enterobacter_spp')
            elif 'klebsiella_pneumoniae' not in self.pathogen_intercept_winner and self.syndrome_no_claims == 'meningitis':
                self.pathogen_intercept_winner.append('klebsiella_pneumoniae')
            elif 'enterobacter_spp' not in self.pathogen_intercept_winner and self.syndrome_no_claims == 'lower_respiratory_infection':
                self.pathogen_intercept_winner.append('enterobacter_spp')
            elif 'mycobacterium_others' in model and self.syndrome_no_claims == 'skin_infection':
                self.pathogen_intercept_winner.remove(model)
            elif 'chlamydia_spp' in model and self.syndrome_no_claims == 'lower_respiratory_infection' and self.model_name == 'community':
                self.pathogen_intercept_winner.remove(model)
            elif 'staphylococcus_aureus' in model and self.syndrome_no_claims == 'lower_respiratory_infection':
                self.pathogen_intercept_winner.remove(model)
            elif 'streptococcus_pneumoniae' in model and self.syndrome_no_claims == 'lower_respiratory_infection':
                self.pathogen_intercept_winner.remove(model)
            elif 'acinetobacter_baumannii' in model and self.syndrome_no_claims == 'lower_respiratory_infection' and self.model_name == 'community':
                self.pathogen_intercept_winner.remove(model)
            elif 'acinetobacter_baumannii' in model and self.syndrome_no_claims == 'lower_respiratory_infection' and self.model_name == 'hospital':
                self.pathogen_intercept_winner.remove(model)
        
        self.log.info(f"These will be in intercept instead of individual:{self.pathogen_intercept_winner}")

        for m in ['individual', 'intercept', 'family', 'other']:
            try:
                self.box_plot_by_pathogen_ordered(model=m)
                self.log.info(f"- - Box plot for {m} MODEL generated!")
            except:
                self.log.warning(f"- - Failed to generate box plot for {m} MODEL!")
        
        filename = "FILEPATH"
        with open(filename, "w") as fp:
            json.dump(self.save_syndrome, fp, default=str, indent=4, sort_keys=True)
        self.log.info(f"- - Saved coefficients for models at {filename}\n")
        
        # GENERATE PREDICTIONS
        preds = self.r_plots_data()
        self.final_df = preds
        
        # PATHOGEN RANK TABLES
        try:
            self.log.info(f"\tMaking pathogen rank table,")
            self.make_pathogen_rank_table(preds)
        except Exception as e:
            self.log.critical(f"{e}")
        
        # CREATING INTERACTIVE GRAPHS FOR VETTING
        self.log.info(f"Innitiating interactive graph for individual pathogens")
        if self.limited:
            pass
        else:
            saveloc = self.create_individual_interactive_graphs(preds, allmodel=False)
            self.log.info(f"- ALL INDIVIDUAL graph saved to: {saveloc}\n")
            self.create_interactive_pathogen_graphs(preds)
            self.log.info(f"- ALL graphs saved to: {saveloc}\n")

        saveloc = self.create_interactive_age_graphs(preds)
        self.log.info(f"- AGE graphs saved to: {saveloc}\n")
        
        
        self.log.info("GENERATING OUTLIERS")
        self.filter_outliers()
        self.log.info("GENERATING RESULTS")
        save_path = "FILEPATH"
        self.result_df = self.gen_results(limited=self.limited)
        self.result_df.to_csv(save_path, index=False)

        self.log.info(f"\t\t **** FINISHED!!!")

    def filter_outliers(self, fraction=25):
        """
        """
        def get_percent(df, fraction):
            threshold = np.percentile(df.standardized_res , [fraction,100-fraction])
            perc_df = df.loc[(df.standardized_res <= threshold[0])|(df.standardized_res >= threshold[1])]
            aggregations = {
                'haqi': 'mean',
                'standardized_res': 'mean',
                'cases': 'sum',
                'deaths': 'sum'
                }
            cases_sum = perc_df.groupby(['spec_source', 'pathogen', 'agg_age_id', 'year_id']).agg(aggregations).reset_index()
            cases_sum.rename(columns={
                'haqi': 'mean_haqi',
                'standardized_res': 'mean_standardized_res',
                'cases': 'outlier_sum_cases',
                'deaths': 'outlier_sum_deaths'
                }, inplace=True)
            return cases_sum
        self.outliers = {}
        for model in [k for k in self.output.keys() if k in ['individual', 'intercept', 'family']]:
            dataframe = self.output[model]
            totals_df = dataframe.groupby(['spec_source', 'pathogen', 'agg_age_id', 'year_id']).agg({'cases': 'sum',
                                                                            'deaths': 'sum'
                                                                        }).reset_index()
            totals_df.rename(columns={
                'cases': 'total_cases',
                'deaths': 'total_deaths'}, inplace=True)
            if model in ['intercept', 'family', 'other']:
                save_df = get_percent(dataframe, fraction)
            else:
                outlier_list = []
                for pat, age in dataframe[['pathogen', 'agg_age_id']].drop_duplicates().itertuples(index=False):
                    tempdf = dataframe.loc[(dataframe.pathogen==pat)&(dataframe.agg_age_id==age)]
                    self.log.info(f"Working on outliers for {pat} {age} for {model}.")
                    top_10_percent_df = get_percent(tempdf, fraction)
                    self.log.info(f"{top_10_percent_df.columns}")
                    outlier_list.append(top_10_percent_df)
                save_df = pd.concat(outlier_list, ignore_index=True)
                self.log.info(f"After individual for loop and concat {save_df.columns}")
            save_df = pd.merge(save_df,
                            totals_df[['spec_source', 'pathogen', 'agg_age_id', 'year_id', 'total_cases', 'total_deaths']],
                            on=['spec_source', 'pathogen', 'agg_age_id', 'year_id'],
                            how='left')
            save_df['outlier_cfr'] = save_df['outlier_sum_deaths']/save_df['outlier_sum_cases']
            save_df['observed_cfr'] = save_df['total_deaths']/save_df['total_cases']
            save_df = save_df.sort_values(by='mean_standardized_res', ascending=False)
            file_path = "FILEPATH"
            save_df.to_csv(file_path, index=False)
            self.log.info(f"Saved top {fraction}% outliers at {file_path}.")
            self.log.info(f"Available with self.outliers[model], with columns: {save_df.columns}")
            self.outliers[model] = save_df
        self.log.info(f"DONE WITH OUTLIERS {self.syndrome}.")
        return None

    def get_vacc_model(self, data_with_vacc, vacc_pat):
        univariate_lr = LinearRegression(fit_intercept=False)
        if vacc_pat == 'streptococcus_pneumoniae':
            coverage_col = 'PCV3_coverage_prop'
        elif vacc_pat == 'haemophilus_influenzae':
            coverage_col = 'Hib3_coverage_prop'
        pcvdata = data_with_vacc.loc[data_with_vacc.pathogen == vacc_pat]
        x = pcvdata.haqi.to_numpy().reshape(-1, 1)
        univariate_lr.fit(x, pcvdata[coverage_col].to_numpy())
        self.vacc_preds[vacc_pat] = {'slope': univariate_lr.coef_[0],
                                                    'r_value': univariate_lr.score(x, pcvdata[coverage_col].to_numpy()),
                                                    'pcv_lin_model': univariate_lr,
                                                    }
        self.log.info(f"Vaccine model for {vacc_pat} saved to self.vacc_preds[{vacc_pat}].")
        return None

    def make_point_preds(self, preds, limited=False):
        swapped_dict = {v: k for k, v in self.age_names_dict.items()}
        preds['agg_age_id'] = preds.agg_age_id.replace(swapped_dict)
        pat_list = self.estimate_pathogen_dict[self.syndrome_no_claims]
        age_col = self.model_ages
        lh = get_current_location_hierarchy(location_set_version_id=self.conf.get_id('location_set_version')).query("level == 3")
        locs = lh.location_id.unique().tolist()
        if limited:
            years = [1990, 2000, 2010 , 2021]
        else:
            years = np.arange(1980, 2023, 1).tolist()
        index = pd.MultiIndex.from_product(
                    [locs, age_col, years],
                    names=['location_id', 'agg_age_id', 'year_id',])
        square_df = pd.DataFrame(index=index).reset_index()
        square_df = square_df.assign(sex_id=3)
        square_df = merge_covariate(square_df, 'haqi', release_id=self.conf.get_id("release"),)
        square_df['haqi'] = square_df.haqi.round(3)
        finaldf = pd.DataFrame()
        square_df['haqi'] = square_df['haqi'].astype(float)
        square_df['agg_age_id'] = square_df['agg_age_id'].astype(int)
        preds['haqi'] = preds['haqi'].astype(float)
        preds['agg_age_id'] = preds['agg_age_id'].astype(int)
        for pathogen in preds.pathogen.unique():
            if pathogen in ['all_fungi', 'all_virus', 'all_bacteria', 'all_parasite']:
                preds = preds.loc[preds.pathogen != pathogen]
            else:
                tempdf = preds.loc[preds.pathogen == pathogen].copy()
                unique_models = tempdf.model.unique()
                if len(unique_models) > 1 and 'intercept' in unique_models:
                    excluded_index = tempdf.loc[(tempdf.pathogen == pathogen) & (tempdf.model == 'intercept')].index
                    tempdf = tempdf.drop(index = excluded_index)
                if pathogen in pat_list:
                    finaldf = pd.concat([finaldf,pd.merge(square_df, tempdf,
                                                          on=['haqi','agg_age_id'])], ignore_index=True)
        not_predicted_yet = list(set(pat_list)-set(finaldf.pathogen.unique()))
        self.log.info(f"Pathogens not predicted yet - USING OTHER MODEL: {not_predicted_yet}")
        for pat in not_predicted_yet:
            tempdf = preds.loc[preds.pathogen == 'all'].copy()
            tempdf['pathogen'] = pat
            finaldf = pd.concat([finaldf,pd.merge(square_df, tempdf,
                                                  on=['haqi','agg_age_id'])], ignore_index=True)
        dropcols = [c for c in ['ICU_ICU_only', 'intercept', 'model', 'PCV3_coverage_prop', 'Hib3_coverage_prop'] if c in finaldf.columns]
        finaldf = finaldf.drop(columns=dropcols)
        return finaldf
    
    def gen_results(self, limited=True, m_type=False):
        def process_binom_model(model, sq_df, pathogenlist):
            model_columns = list(self.binom_df[model].columns)
            template_df = sq_df.copy()
            if any([c.startswith('ICU_') for c in model_columns]):
                template_df['ICU_ICU_only'] = 0
            if any([c.startswith('hosp_') for c in model_columns]):
                template_df['hosp_'+self.get_value("hosp")] = 1
            def expand_list_to_cols(df, col_list, prefix):
                temp = df.copy()
                temp[prefix+'_'+str(self.save_syndrome[model+'_binom']['references'][prefix])] = 1
                if prefix == 'pathogen':
                    col_list = [prefix+'_'+col for col in col_list]
                for col in col_list:
                    mergedf = df.copy()
                    mergedf[col] = 1
                    temp = pd.concat([temp, mergedf], ignore_index=True)
                return temp.fillna(0)
            ages_in_model = [col for col in model_columns if 'agg_age_id' in col]
            template_df = expand_list_to_cols(template_df, ages_in_model, 'agg_age_id')
            model_pathogens = [str(col).replace('pathogen_', '') for col in model_columns if 'pathogen_' in col]
            template_df = expand_list_to_cols(template_df, model_pathogens, 'pathogen')
            model_pathogens = model_pathogens + [self.save_syndrome[model+'_binom']['references']['pathogen']]
            if model not in ['intercept', 'family', 'other']:
                model_pathogens = [model]
            pathogens_not_estimated = [pat for pat in model_pathogens if pat not in pathogenlist]
            estimated_pathogens = [pat for pat in model_pathogens if pat in pathogenlist]
            self.log.info(f"\t*Pathogens predicted by {model} model: {model_pathogens}")
            template_df['intercept'] = 1
            # self.log.info(f"\t*Template dataframe updated for with columns {template_df.columns} model, and shape {template_df.shape}.")
            checklist = estimated_pathogens+pathogens_not_estimated
            template_df = self.merge_vacc_data(template_df, checklist)
            cfr_results = self.regmod_binom_model[model].predict(template_df)
            not_included = ['pathogen_'+p for p in pathogens_not_estimated]
            try:
                mask = (cfr_results[not_included] == 1).any(axis=1)
                filtered_cfr_results = cfr_results[~mask]
                filtered_cfr_results = filtered_cfr_results.drop(columns=not_included)
            except:
                self.log.info(f"Could not filter out {not_included} from {cfr_results.columns} for -{model}-.")
                filtered_cfr_results = cfr_results
            return filtered_cfr_results
        lh = get_current_location_hierarchy(
            location_set_version_id= self.conf.get_id('location_set_version')).query("level == 3")
        locs = lh.location_id.unique().tolist()
        if limited:
            years = [1990, 2005, 2021]
        else:
            years = np.arange(1980, 2023, 1).tolist()
        all_dfs = []
        index = pd.MultiIndex.from_product(
            [locs, years],
            names=['location_id', 'year_id',])
        square_df = pd.DataFrame(index=index).reset_index()
        square_df = merge_covariate(square_df, 'haqi', release_id=self.conf.get_id("release"),)
        
        if self.syndrome_no_claims not in self.estimate_pathogen_dict.keys():
            self.log.critical(f"Syndrome {self.syndrome} not in pathogens to assess by syndrome file, and has no pathogen list.")
            return None
        pathogenlist = self.estimate_pathogen_dict[self.syndrome_no_claims]+['all']
        individual_pats = []
        if m_type:
            run_list = [m_type]
        else:
            run_list = self.binom_df.keys()
        for model in run_list:
            if model not in ['intercept', 'family', 'other']:
                if model in self.pathogen_intercept_winner:
                    continue
                indv_df = process_binom_model(model, square_df, pathogenlist)
                self.log.info(f"\t*Processing:{model}")
                indv_df = self.collapse_col('intercept_binom', indv_df,'agg_age_id')
                indv_df['pathogen'] = model
                self.log.info(f"\t*Processed:{indv_df.columns}\t{indv_df.agg_age_id.dtype}")
                indv_df['agg_age_id'] = indv_df['agg_age_id'].astype(int)
                try:
                    self.log.info(f"\t*{max(indv_df.loc[indv_df.agg_age_id == 239].p)}.")
                    tempval = abs(max(indv_df.loc[indv_df.agg_age_id == 239].p)- min(indv_df.loc[indv_df.agg_age_id == 239].p))
                    self.log.info(f"\t*{model} model has max-min of {tempval}\n")
                    if indv_df.p.max() > 0.7:
                        self.log.warning(f"**** {model} model has predictions > 0.8 CFR")
                    elif tempval < 0.01:
                        self.log.warning(f"**** {model} model IS FLAT with absolute difference of {indv_df.p.max()-indv_df.p.min()}")
                    else:
                        all_dfs.append(indv_df)
                        individual_pats.append(model)
                except Exception as e:
                    self.log.critical(f"Error in {model} model: {e}")
                    self.log.critical(f"\t* ages in model: {indv_df.agg_age_id.unique()}")
                    self.log.critical(f"\t* possible issues with prediction {indv_df.p.describe()}")
                    continue
        if len(individual_pats) > 0:
            self.log.info(f"\t*Processed:\n{individual_pats}\n\t***Pathogens as individual.\n")
        else:
            self.log.critical(f"\t* There are no individual models to produce predictions.") 
        if 'intercept' in run_list:
            int_df = process_binom_model('intercept', square_df, pathogenlist)
            int_df = self.collapse_col('intercept_binom', int_df,'agg_age_id')
            int_df = self.collapse_col('intercept_binom', int_df,'pathogen')
            self.log.info(f"\t*Processed:\n{int_df.columns}\n\t\tAS PREDICTIONS for intercept model.")
            if any(int_df.p > 0.7):
                high_cfr_pats = int_df.loc[int_df.p > 0.7].pathogen.unique()
                self.log.warning(f"**** {high_cfr_pats} model has predictions > 0.7 CFR")
                for pat in high_cfr_pats:
                    int_df = int_df.loc[int_df.pathogen != pat]
            else:
                self.log.info(f"\t* No predictions > 0.7 for intercept model.")
            intercept_pathogens = int_df.pathogen.unique()
            self.log.info(f"\t*Pathogens found in intercept model:\n{intercept_pathogens}\n")
            intercept_pats = list(set(intercept_pathogens)-set(individual_pats))
            self.log.info(f"\t*Pathogens found only in intercept model:\n{intercept_pats}\n")
            for int_pat in intercept_pats:
                tempdf = int_df.loc[int_df.pathogen == int_pat].copy()
                all_dfs.append(tempdf)
            remaining_pats = [pat for pat in pathogenlist if pat not in (intercept_pats+individual_pats)]
            self.log.info(f"\t*Processed:\n{intercept_pats}\n\tPathogens from intercept model.")
        else:
            remaining_pats = [p for p in pathogenlist if p not in individual_pats]
            self.log.critical(f"Could not process pathogens from pathogen list from intercept model. Pathogen list:\n{pathogenlist}")
        # ASSIGN TO FAMILY IF POSSIBLE
        try:
            if 'family' in run_list:
                in_family_model = [col.replace('pathogen_','') for col in self.binom_df['family'].columns if 'pathogen_' in col
                                ] + [self.reference_dict['family']['pathogen']]
                fam_preds = process_binom_model('family', square_df, in_family_model)
                fam_preds = self.collapse_col('family_binom', fam_preds,'agg_age_id')
                fam_preds = self.collapse_col('family_binom', fam_preds,'pathogen')
                iterateonthis = remaining_pats.copy()
                for pat in iterateonthis:
                    if pat not in self.pat_family_dict.keys():
                        self.log.warning(f"\t*Pathogen:{pat} not in pathogen to family dictionary, will use other/all model.")
                        continue
                    fam = self.pat_family_dict[pat]
                    if fam in in_family_model:
                        tempdf = fam_preds.loc[fam_preds.pathogen == fam].copy()
                        tempdf['pathogen'] = pat
                        all_dfs.append(tempdf)
                        remaining_pats.remove(pat)
                    else:
                        self.log.warning(f"Pathogen:{pat} with family:{fam} not in family model.")
        except Exception as e:
            self.log.warning(f"Could not process pathogens from family model. Error: {e}")
        if len(remaining_pats) > 0:
            # ASSIGN ALL UNMODELED PATHOGENS TO 'OTHER' MODEL
            self.log.info(f"\t* Processing {remaining_pats}\n\t\t* As 'other' binomial model.")
            all_preds = process_binom_model('other', square_df, ['all'])
            all_preds = self.collapse_col('other_binom', all_preds,'agg_age_id')
            for pat in remaining_pats:
                tempdf = all_preds.copy()
                tempdf['pathogen'] = pat
                all_dfs.append(tempdf)
        final_df = pd.concat(all_dfs, ignore_index=True)
        patnotinlist = [p for p in final_df.pathogen.unique() if p not in pathogenlist]
        if len(patnotinlist) > 0:
            self.log.warning(f"PATHOGEN NOT IN LIST:\n\t\t{patnotinlist}.")
        if len([p for p in pathogenlist if p not in final_df.pathogen.unique()]) > 0:
            self.log.critical(f"\t**** NO CFR PREDICTION FOR {[p for p in pathogenlist if p not in final_df.pathogen.unique()]}.")
        final_df.rename(columns={'p': 'predict', 'agg_age_id':'age_group_id'}, inplace=True)
        if self.get_value("hosp"):
            final_df = final_df.assign(hosp=self.get_value("hosp"))
        final_df = final_df.assign(sex_id=3)
        keepcols = [c for c in final_df.columns if c in ['location_id', 'year_id', 'haqi', 'age_group_id','sex_id', 'hosp','pathogen', 'predict']]
        return final_df[keepcols]

    def make_pathogen_rank_table(self, pred_df):
        excel_file = "FILEPATH"
        writer = pd.ExcelWriter(excel_file, engine='xlsxwriter')
        df = pred_df.drop(columns=['haqi', 'intercept'], axis=1)
        for m_type in df.model.unique():
            tempdf = df.loc[df.model == m_type].copy()
            if m_type == 'individual':
                moved_individual = tempdf.loc[tempdf.pathogen.isin(self.pathogen_intercept_winner)]
                pivot_moved = moved_individual.pivot_table(index='pathogen', columns='agg_age_id', values='predict', aggfunc='mean').reset_index()
                pivot_moved.to_excel(writer, sheet_name='moved_individual', index=False)
                tempdf = tempdf.loc[~tempdf.pathogen.isin(self.pathogen_intercept_winner)]
            pivot_df = tempdf.pivot_table(index='pathogen', columns='agg_age_id', values='predict', aggfunc='mean').reset_index()
            if m_type in ['intercept', 'individual']:
                removed = [p for p in pivot_df.pathogen.unique() if p not in self.estimate_pathogen_dict[self.syndrome_no_claims]]
                missing = [p for p in self.estimate_pathogen_dict[self.syndrome_no_claims] if p not in pivot_df.pathogen.unique()]
                self.log.warning(f"\t* Pathogens NOT ESTIMATED:\n\t{removed}")
                self.log.warning(f"\t* Pathogens moved to either FAMILY or OTHER/ALL: {missing}")
                pivot_df = pivot_df.loc[pivot_df.pathogen.isin(self.estimate_pathogen_dict[self.syndrome_no_claims])]
            numeric_cols = pivot_df.select_dtypes(include=['number'])
            pivot_df['average_cfr'] = numeric_cols.mean(axis=1)
            pivot_df = pivot_df.sort_values(by='average_cfr', ascending=False)
            pivot_df.to_excel(writer, sheet_name=m_type, index=False)
        writer.save()
        self.log.info(f"Pathogen rank table saved to {excel_file}.")

    def generate_scores_table(self, pat):
        def calculate_metrics_by_group(df, group_by_column='agg_age_id'):
            # This function calculates RMSE and R² for each group specified by `group_by_column`
            scores_by_group = {}
            for group, group_df in df.groupby(group_by_column):
                rmse = np.sqrt(mean_squared_error(group_df['observed'], group_df['predict']))
                r_squared = r2_score(group_df['observed'], group_df['predict'])
                scores_by_group[group] = {'RMSE': rmse, 'R2': r_squared}
            return scores_by_group
        def calculate_overall_scores(df):
            # This function calculates the overall RMSE and R² for the dataset
            overall_rmse = np.sqrt(mean_squared_error(df['observed'], df['predict']))
            overall_r_squared = r2_score(df['observed'], df['predict'])
            return {'Overall RMSE': overall_rmse, 'Overall R2': overall_r_squared}
        df1 = self.output[pat]
        intercept_df = self.output['intercept']
        df2 = intercept_df.loc[intercept_df.pathogen == pat]
        df1_scores_by_group = calculate_metrics_by_group(df1)
        df2_scores_by_group = calculate_metrics_by_group(df2)
        df1_overall_scores = calculate_overall_scores(df1)
        df2_overall_scores = calculate_overall_scores(df2)
        combined_scores = {}
        for age_id in set(list(df1_scores_by_group.keys()) + list(df2_scores_by_group.keys())):
            combined_scores[age_id] = {
                'df1_RMSE': df1_scores_by_group.get(age_id, {}).get('RMSE', None),
                'df1_R2': df1_scores_by_group.get(age_id, {}).get('R2', None),
                'df2_RMSE': df2_scores_by_group.get(age_id, {}).get('RMSE', None),
                'df2_R2': df2_scores_by_group.get(age_id, {}).get('R2', None)
            }
        combined_scores['Overall'] = {
            'df1_RMSE': df1_overall_scores.get('Overall RMSE', None),
            'df1_R2': df1_overall_scores.get('Overall R2', None),
            'df2_RMSE': df2_overall_scores.get('Overall RMSE', None),
            'df2_R2': df2_overall_scores.get('Overall R2', None)
        }
        scores_df = pd.DataFrame.from_dict(combined_scores, orient='index').reset_index().rename(columns={'index': 'agg_age_id'})
        scores_df = scores_df[['agg_age_id', 'df1_RMSE','df2_RMSE', 'df1_R2', 'df2_R2']]
        scores_df['by_RSME'] = np.where(scores_df.df1_RMSE<scores_df.df2_RMSE, 'individual', 'intercept')
        scores_df['by_R2'] = np.where(scores_df.df1_R2>scores_df.df2_R2, 'individual', 'intercept')
        return scores_df

    def pick_model(self):
        intercept_winner = []
        scores_all = []
        for pat in [p for p in self.output.keys() if p not in ['other','family','intercept','individual']]:
            scores_table = self.generate_scores_table(pat)
            if len(scores_table) % 2 != 0:
                winner = scores_table['by_R2'].value_counts().idxmax()
            else:
                value_counts = scores_table['by_R2'].value_counts()
                if len(value_counts) == 2 and value_counts.iloc[0] == value_counts.iloc[1]:
                    winner = scores_table.loc[scores_table['agg_age_id'] == 'Overall', 'by_R2'].iloc[0]
                else:
                    winner = scores_table['by_R2'].value_counts().idxmax()
            if winner == 'intercept':
                intercept_winner.append(pat)
            self.log.info(f"===={pat}====\n{scores_table}\nWINNER IS:{winner}")
            scores_table['pathogen'] = pat
            scores_all.append(scores_table)
        save_scores = pd.concat(scores_all, ignore_index=True)
        save_scores.to_csv("FILEPATH", index=False)
        return intercept_winner
