import pandas as pd
import numpy as np
import argparse
import re
import sys
from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import get_ages, get_map_version, pretty_print
from cod_prep.utils import (report_if_merge_fail, print_log_message,
                            clean_icd_codes)
from cod_prep.utils.formatting import ages
from mcod_prep.mcod_mapping import MCoDMapper
from amr_prep.utils.amr_io import save_amr_data
from amr_prep.utils import assign_breakpoints

CONF = Configurator()
SOURCE = "SOURCE"
L_DIR = "FILEPATH"
repo_dir = "FILEPATH"
MAP_DIR = "FILEPATH"

breakpt_dir = "FILEPATH"
intermediate_dir = "FILEPATH"
raw_dir = "FILEPATH"

def create_age_group_id(df):

    df.loc[df.age < 0, "age"] = np.NaN

    df['age_unit'] = 'year'
    df.rename(columns={'Age': 'age'}, inplace=True)
    age_formatter = ages.PointAgeFormatter()
    df = age_formatter.run(df)

    assert df['age_group_id'].notnull().values.all()

    return df


def create_sex_id(df):
    df['sex_id'] = df['sex'].apply(
        lambda x: 1 if x == 'male' else (2 if x == 'female' else 9)
    )
    print(
        f"{(len(df.loc[df['sex_id'] == 9]) / len(df)) * 100}% of rows are missing sex")
    report_if_merge_fail(df, 'sex_id', 'sex')
    return df


def mark_admissions_deaths(df):
    discharge_map = pd.read_csv("FILEPATH")
    discharge_dict = discharge_map.set_index("discharge_disp")[
        "Died"].to_dict()

    df["admissions"] = 1

    df["deaths"] = df["discharge_disp"].map(discharge_dict)
    df.loc[(df.deaths.isnull()) & (df.discharge_disp.str.contains(
        "hospice", flags=re.IGNORECASE)), "deaths"] = 1
    df.loc[df.discharge_disp.isin(["nan", "UNKNOWN", "-"]), "deaths"] = np.nan
    return df


def apply_pathogen_severity_map(df, named_cols):
    pathogen_severity = pd.read_excel("FILEPATH")
    pathogen_dict = pathogen_severity.set_index("pathogen_name")["rank"].to_dict()

    df[named_cols] = df[named_cols].apply(lambda x: x.map(
        pathogen_dict)).fillna(100).astype(int)

    df["severe_pathogen_int"] = df[named_cols].apply(lambda x: x.min(), axis=1)

    reverse_dict = dict([value, key] for key, value in pathogen_dict.items())
    df["severe_pathogen"] = df["severe_pathogen_int"].map(
        reverse_dict)

    df[named_cols] = df[named_cols].apply(lambda x: x.map(
        reverse_dict))

    df["severe_pathogen_int"] = df["severe_pathogen_int"].replace(100, np.NaN)

    return df


def fix_multiple_deaths(df):
    date_map = pd.read_csv("FILEPATH")
    df = mark_admissions_deaths(df)

    date_map.rename(columns={"date_id": "discharge_date_id",
                             "date": "discharge_date"}, inplace=True)
    df = df.merge(date_map, on="discharge_date_id", how="left")
    df["discharge_date"] = pd.to_datetime(df.discharge_date)
    report_if_merge_fail(df, "discharge_date_id", "date")
    null = df.loc[df['patient_id'].isna(), ]
    notnull = df.loc[~df['patient_id'].isna(), ]
    notnull["max"] = notnull.groupby("patient_id")["discharge_date"].transform("max")
    notnull["deaths"] = np.where(notnull["discharge_date"] == notnull["max"], notnull["deaths"], 0)
    df = pd.concat([null, notnull])

    return df


def map_specimen_source(df):
    specimen_detail = pd.read_excel("FILEPATH")
    specimen_detail.rename(columns={"specimen_source_group": "specimen_source"}, inplace=True)

    df = df.merge(specimen_detail, on=["site_id", "specimen_id", "specimen_source"], how="left")
    df.loc[df.specimen_source_detail.isnull(), "specimen_source_detail"] = df["specimen_source"]
    report_if_merge_fail(df, "specimen_source_detail", "specimen_source")

    df.rename(columns={'specimen_source_detail': 'raw_specimen'}, inplace=True)

    return df


def apply_syndrome_severity_map(df):
    priority_dict = {"peritonitis": 1,
                     "lri": 2,
                     "uti": 3,
                     "sepsis/bacteremia": 4,
                     "cellulitis": 5,
                     "other": 6}

    df['syndrome_num'] = df.groupby('admission_id')['specimen_syndrome'].transform(
        lambda syndromes: [
            f"specimen_syndrome_{i}" for i in range(len(syndromes))]
    )
    syndrome_df = df.pivot(
        index='admission_id', columns='syndrome_num', values='specimen_syndrome'
    ).reset_index()

    syndrome_cols = [x for x in list(syndrome_df) if "specimen_syndrome" in x]
    syndrome_df[syndrome_cols] = syndrome_df[syndrome_cols].apply(lambda x: x.map(
        priority_dict)).fillna(100).astype(int)
    syndrome_df = syndrome_df.assign(severe_syndrome_int=lambda x: x.min(axis=1))

    reverse_dict = dict([value, key] for key, value in priority_dict.items())
    syndrome_df["severe_syndrome"] = syndrome_df["severe_syndrome_int"].map(
        reverse_dict)

    df = df.merge(
        syndrome_df[["admission_id", "severe_syndrome"]], on="admission_id", how="left")
    drop_synds = df[df["specimen_syndrome"] != df["severe_syndrome"]].index
    df.drop(drop_synds, inplace=True)

    df.drop(columns="severe_syndrome", axis=1, inplace=True)

    return df


def create_pathogen_cols(df):

    df.rename(columns={'org_name': 'raw_pathogen'}, inplace=True)
    df['raw_pathogen'] = df['raw_pathogen'].str.strip()

    return df


def mark_estimated_pathogens(df, cols):
    df["pathogen"] = np.where(df[cols].isnull().all(axis=1), "none", "other")
    df["pathogen"] = np.where(df["severe_pathogen"].notnull(),
                              df["severe_pathogen"], df["pathogen"])

    return df


def map_cause_cols(df):
    print_log_message("creating pathogens cols")
    df = create_pathogen_cols(df)
    diag_df = prep_diagnosis()

    df = df.merge(diag_df, on="admission_id", how="left")
    df.rename(columns={"primary_diagnosis_code": "cause"}, inplace=True)

    print_log_message("mapping specimen source col")
    df = map_specimen_source(df)

    print_log_message("fixings patients with multiple 'deaths'")
    df = fix_multiple_deaths(df)

    for col in [x for x in list(df) if "cause" in x]:
        df[col] = df[col].fillna("none")
    df['cause'] = clean_icd_codes(df['cause'], remove_decimal=True)
    df.loc[~df.cause.str.match(r"[A-Z]\d{2,5}"), 'cause'] = "none"


    return df


def prep_diagnosis():
    df = pd.read_csv(
        "FILEPATH",
        usecols=['admission_id', 'primary_diagnosis', 'diagnosis_code'],
        dtype={'diagnosis_code': object}
    )
    df['diagnosis_code'] = clean_icd_codes(df['diagnosis_code'], remove_decimal=True)
    df.loc[~df.diagnosis_code.str.match(r"[A-Z]\d{2,5}"), 'diagnosis_code'] = "none"
    df['cause_num'] = df.groupby('admission_id')['diagnosis_code'].transform(
        lambda codes: [f"multiple_cause_{i}" for i in range(len(codes))]
    )
    df.loc[df['primary_diagnosis'], 'cause_num'] = 'cause'
    df = df.pivot(
        index='admission_id', columns='cause_num', values='diagnosis_code'
    ).reset_index().fillna('none')
    df.drop(columns="cause", inplace=True)
    return df


def map_hosp(df):
    date_map = pd.read_csv("FILEPATH")

    admit_date_map = date_map.rename(columns={"date_id": "admit_date_id",
                                              "date": "admit_date"})
    collected_date_map = date_map.rename(columns={"date_id": "collected_date_id",
                                                  "date": "collected_date"})
    dis_date_map = date_map.rename(columns={"date_id": "discharge_date_id",
                                              "date": "dis_date"})

    df = df.merge(admit_date_map, on='admit_date_id', how='left')
    df = df.merge(collected_date_map, on='collected_date_id', how='left')
    df = df.merge(dis_date_map, on='discharge_date_id', how='left')
    df['date_diff'] = (pd.to_datetime(df['collected_date']) -
                       pd.to_datetime(df['admit_date'])).dt.days

    df['hosp'] = 'hospital'
    df.loc[df['date_diff'] <= 2, 'hosp'] = 'community'
    df.loc[df['date_diff'].isna(), 'hosp'] = 'unknown'
    df['days_prior'] = (pd.to_datetime(df['collected_date']) -
                      pd.to_datetime(df['admit_date'])).dt.days

    df['days_infection'] = (pd.to_datetime(df['dis_date']) -
                      pd.to_datetime(df['admit_date'])).dt.days
    return df


def get_ast_values(df):
    
    def extract_number(numeric_str):
    
        equal_pat = r'^(={0,1}[0-9]{1,}\.{0,1}[0-9]{0,})|(={0,1}[0-9]{0,}\.{0,1}[0-9]{1,})'

        less_pat = r'^((<)[0-9]{1,}\.{0,1}[0-9]{0,})|((<)[0-9]{0,}\.{0,1}[0-9]{1,})'
        more_pat =  r'^((>)[0-9]{1,}\.{0,1}[0-9]{0,})|((>)[0-9]{0,}\.{0,1}[0-9]{1,})'

        less_equal = r'^((<=)[0-9]{1,}\.{0,1}[0-9]{0,})|((<=)[0-9]{0,}\.{0,1}[0-9]{1,})'
        more_equal = r'^((>=)[0-9]{1,}\.{0,1}[0-9]{0,})|((>=)[0-9]{0,}\.{0,1}[0-9]{1,})'



        num_slice = r'([0-9]{1,}\.{0,1}[0-9]{0,})|([0-9]{0,}\.{0,1}[0-9]{1,})'
        divide_pat = r'^(<|>|(>=)|(<=)|=){0,1}[0-9]{0,5}\.{0,1}[0-9]{0,3}\/{1}[0-9]{0,5}\.{0,1}[0-9]{0,3}$'

        numeric_str = str(numeric_str)

        if re.match(divide_pat, numeric_str) and (re.search(num_slice, numeric_str)):
            number_1 = numeric_str.split('/')[0]
            try: 
                number_2 = float(numeric_str.split('/')[1])
            except:
                number_2 = 1

            try:
                number_1 = float(number_1[re.search(num_slice, number_1).start():])
            except:
                number_1 =np.nan
            number = number_1
        
            if re.search(r'<=', numeric_str):
                number_val = 'less_equal'

            elif re.search(r'>=', numeric_str):
                number_val = 'more_equal'

            elif re.search(r'>', numeric_str):
                number_val = 'more'

            elif re.search(r'<', numeric_str):
                number_val = 'less'
            else:
                number_val = 'equal'

        elif re.match(equal_pat, numeric_str) and (re.search(num_slice, numeric_str)):  
            number =numeric_str[re.search(num_slice, numeric_str).start():re.search(num_slice, numeric_str).end()]      
            number = float(number)
            number_val = 'equal'

        elif re.match(less_pat, numeric_str) and re.search(num_slice, numeric_str):

            number = float(numeric_str[re.search(num_slice, numeric_str).start():re.search(num_slice, numeric_str).end()])
            number_val = 'less'

        elif re.match(more_equal, numeric_str) and (re.search(num_slice, numeric_str)):

            number = float(numeric_str[re.search(num_slice, numeric_str).start():re.search(num_slice, numeric_str).end()]) 
            number_val = 'more_equal'

        elif re.match(less_equal, numeric_str) and re.search(num_slice, numeric_str):

            number = float(numeric_str[re.search(num_slice, numeric_str).start():re.search(num_slice, numeric_str).end()])
            number_val = 'less_equal'

        elif re.match(more_pat, numeric_str) and (re.search(num_slice, numeric_str)):

            number = float(numeric_str[re.search(num_slice, numeric_str).start():re.search(num_slice, numeric_str).end()]) 
            number_val = 'more'

        else:
            number = np.nan
            number_val = np.nan


        return number, number_val             


    df['value'] = df.raw_result.str.strip()
    df['tmp'] = df.value.apply(extract_number)
    df['number'] = df['tmp'].apply(lambda x: x[0])
    df['number_val'] = df['tmp'].apply(lambda x: x[1])

    df['raw_method'] = df.method

    df['method'] = np.nan

    mic_methods = ['MIC','ETEST', 'MIC-Systemic', 'MIC-Urine','Etest','BACTEC','ESBL','AFB MIC','AFB','DIL', 'AFB SUSC','BLAC','Beta Lactamase']

    df.loc[df.raw_method.isin(mic_methods) , 'method']='MIC'

    explore_methods = ['UNK', 'KB', 'BP','OTHER','AFB']

    def decide_mathod(numeric_str):   

        numeric_str = str(numeric_str)
        disk_pat = r'(MM)|(mm)'

        mic_pat = r'(ug/mL)|(MCG/ML)|MIC'


        if re.search(disk_pat, numeric_str):
            return 'DISK'
        elif re.search(mic_pat, numeric_str):
            return 'MIC'
        else:
            return np.nan

    df.loc[(df.raw_method.isin(explore_methods))&(~df.value.isna()), 'method']= df.value.apply(decide_mathod)

    df.loc[df.method.isna()&(df.number<1), 'method']='MIC'
    df.loc[df.method.isna()&(df.number>50), 'method']='MIC'
    
    return df

def add_disk_doses(df):
    
    start_len = len(df)
    
    brm = pd.read_csv("FILEPATH")
    pathogen_map = pd.read_csv("FILEPATH")
    pathogen_map.drop_duplicates('raw_pathogen', inplace=True)

    brm = brm[brm.method=='DISK']
    brm = brm[['pathogen', 'antibiotic','site_category','guideline', 'disk_dose_1']]
    grp = brm.groupby(['pathogen', 'antibiotic','site_category','disk_dose_1']).count().reset_index()


    dup = grp[grp.duplicated(['pathogen','antibiotic', 'site_category'])]

    dup['unique_disk_zone'] =0

    dup = dup.drop(['disk_dose_1','guideline'], axis=1)

    grp = dup.merge(grp, how= 'right', on=['pathogen', 'antibiotic','site_category'])

    safe_combos = grp.loc[grp.unique_disk_zone!=0]

    
    disk_doses = safe_combos[['pathogen','antibiotic','disk_dose_1','site_category']].merge(brm)

    disk_doses = disk_doses[['pathogen','antibiotic','disk_dose_1','site_category']]
    disk_doses = disk_doses.drop_duplicates(['pathogen','antibiotic','site_category'])
    
   
    df_disk = df[df.method=='DISK']   
    df_n_disk = df[df.method!='DISK']

    df_disk = df_disk.merge(pathogen_map, how='left', on='raw_pathogen', validate = 'm:1')

    df_disk= df_disk.merge(disk_doses, how= 'left', left_on=['mapped_pathogen','raw_antibiotic','site_category'], right_on=['pathogen', 'antibiotic','site_category'], indicator=True, validate = "m:1")
    

    
    df_disk_1 = df_disk[df_disk._merge=='both']
    df_disk_2 = df_disk[df_disk._merge=='left_only']
   
    df_disk_2 = df_disk_2.drop(['_merge','disk_dose_1','pathogen','antibiotic'], axis=1)
   
    df_disk_2 = df_disk_2.merge(disk_doses, how= 'left', left_on=['pathogen_name_2','raw_antibiotic','site_category'], right_on=['pathogen', 'antibiotic','site_category'], indicator=True, validate = "m:1")

   
    df_disk = pd.concat([df_disk_1,df_disk_2])
   

    df_disk = df_disk.drop(['mapped_pathogen','pathogen_name_2','pathogen','antibiotic', '_merge','site_category'], axis=1)

    df_disk.rename(columns={'disk_dose_1':'disk_dose'}, inplace=True)

    df= pd.concat([df_disk, df_n_disk])
    assert start_len == len(df)
    
    return df

def more_disk(obs, brk):
    
    if obs >= brk:
        return 'susceptible'
    else:
        return np.nan

def less_disk(obs, brk):
    
    if obs<=brk:
        return 'resistant'
    else:
        return np.nan


def less_equal_disk(obs, brk):
    
    if obs<brk:
        return 'resistant'
    else:
        return np.nan
    
def equal_disk(obs,brk):
    
    if obs >=brk:
        return 'susceptible'
    if obs<brk:
        return 'resistant'
    
def more_mic(obs, brk):
    
    if obs> brk:
        return "resistant"
    else:
        return np.nan

def more_equal_mic(obs,brk):
    
    if obs >= brk:
        return "resistant"
    else:
        return np.nan
    
def less_mic(obs, brk):
    
    if obs <= brk:
        return 'susceptible'
    else:
        return np.nan  

def equal_mic(obs, brk):
    
    if obs <= brk:
        return 'susceptible'
    elif obs > brk:
        return "resistant"

def apply_resistance(s):

    if s.method=='MIC':
        if s.number_val=='equal':
            return equal_mic(s.number, s.breakpoint_S)
        elif s.number_val == 'less' or s.number_val =='less_equal':
            return less_mic(s.number, s.breakpoint_S)
        elif s.number_val == 'more':
            return more_mic(s.number, s.breakpoint_S)
        elif s.number_val == 'more_equal':
            return more_equal_mic(s.number, s.breakpoint_S)
        else:
            return np.nan        
    elif s.method=='DISK':
        if s.number_val=='equal':
            return equal_disk(s.number, s.breakpoint_S)
        elif s.number_val == 'more' or s.number =='more_equal':
            return more_disk(s.number, s.breakpoint_S)
        elif s.number_val == 'less':
            return less_disk(s.number, s.breakpoint_S)
        elif s.number_val == 'less_equal':
            return less_equal_disk(s.number, s.breakpoint_S)
        else:
            return np.nan
    
    else:
        return np.nan

def create_resistance_column_clsi_2023(df):
    
    df['site_category'] = 'other'
    df.loc[df.raw_specimen=='urine', 'site_category'] = 'other'

    df['raw_antibiotic'] = df.raw_antibiotic.str.strip().str.lower()

    csf_rows = ['cefotaxime (csf/mening breakpoints)','ceftriaxone (csf/mening breakpoints)',
                'penicillin / csf', 'benzylpenicillin / csf',
               'benzylpenicillin - csf', 'penicillin - csf']
    df.loc[df.raw_antibiotic.isin(csf_rows), 'site_category'] = 'csf'

    ab_dict = {'rifabutin 1':'rifabutin', 
        'ethionamide 5':'ethionamide',
        'ethambutol 5':'ethambutol',
        'streptomycin 1000':'streptomycin-high',
        'gentamicin 500':'gentamicin',
        'quinupristin/dalfopristin (synercid)':'quinupristin/dalfopristin',
        'cefotaxime (csf/mening breakpoints)':'cefotaxime',
        'doripenem (doribax)':'doripenem',
        'penicillin / csf':'benzylpenicillin',
        'colistin (colistimethate)':'colistin',
        'isoniazid (low level)':'isoniazid',
        'benzylpenicillin e test':'benzylpenicillin',
        'ceftriaxone (csf/mening breakpoints)':'ceftriaxone',
        'cefoxitin screen':'cefoxitin',
        'cefuroxime/axetil (oral)':'cefuroxime axetil',
        'streptomycin high-level':'streptomycin-high',
        'ticarcillin/ca':'ticarcillin/clavulanic acid',
        'cefuroxime/sodium':'cefuroxime',
        'gentamicin high-level': 'gentamicin-high',
        'nitroxolin':'nitroxoline',
        'amphotericin_b':'amphotericin b',
        'tgc':'tigecycline',
        'dorip':'doripenem',
        'methicillin':'meticillin',
        'benzylpenicillin / csf':'benzylpenicillin', 
        'ceftolozane/tazobactam (zerbaxa)':'ceftolozane/tazobactam', 
        'streptomycin':'streptomycin-high', 
        'telavancin (vibativ)':'telavancin', 
        'inh':'isoniazid',   
        'clindamycin d/test':'clindamycin',   
        'pen iv':'benzylpenicillin',  
        'tige':'tigecycline', 
        'penicillin v':'benzylpenicillin',
        'gentamicin 120':'gentamicin',  
        'penicillin / pneum':'benzylpenicillin', 
        'kanamycin 6':'kanamycin',            
        'inh 0.4':'isoniazid', 
        'inh 4':'isoniazid', 
        'inh 0.1':'isoniazid',
        'inh 1':'isoniazid',          
        'ceftazidime/ca':'ceftazidime/avibactam', 
        'ofloxacin 2':'ofloxacin',
        'kanamycin 5':'kanamycin', 
        'cefalothin':'cefalotin', 
        'metronidozole':'metronidazole',
        'rifampin 1':'rifampicin',
        'rifampin':'rifampicin',
        'rifampin 5':'rifampicin', 
        'ethambutol 10':'ethambutol', 
        'polymixin': 'polymyxin b', 
        'rifampin 2':'rifampicin', 
        'ethambutol 7.5':'ethambutol', 
        'mupirocin by pcr':'mupirocin', 
        'tedizolid (sivextro)':'tedizolid', 
        'cefp':'cefprozil',  
        'ethambutol 2.5':'ethambutol',
        'ethambutol 1':'ethambutol', 
        'isoniazid (low level)':'isoniazid', 
        'benzylpenicillin e test':'benzylpenicillin',  
        'sulfamethoxasole':'sulfamethoxazole', 
        'rifa':'rifampicin', 
        'cefaz':'cefazolin',  
        'ethambutol 8':'ethambutol', 
        'strept':'streptomycin-high', 
        'fos':'fosfomycin',
        'penicillin / csf':'benzylpenicillin',
        'cefuroxime/sodium':'cefuroxime',
        'benzylpenicillin - csf':'benzylpenicillin',
        'penicillin - csf':'benzylpenicillin'}

    df['raw_antibiotic'] = df.raw_antibiotic.apply(lambda x: ab_dict[x] if x in ab_dict.keys() else x)
    print(df.raw_antibiotic.unique())

    df = assign_breakpoints.validate_format_antibiotic_names(df)

    df['raw_pathogen'] = df.raw_pathogen.str.replace(' ', '_')
    
    df = add_disk_doses(df)
    

    df = assign_breakpoints.assign_breakpoints(df, 'CLSI 2023')
    df['resistance'] = df.apply(lambda x: apply_resistance(x), axis=1)
    
    df.loc[df.resistance.isna(), 'resistance'] = 'unknown'
    print_log_message(f'{len(df[df.resistance!="unknown"])/len(df)*100} % of the dataset has resistance values assigned')
    
    return df

def clean_SOURCE():
    df = pd.read_csv(
        "FILEPATH",
        encoding='iso-8859-1')
    cols = [col for col in df.columns if 'result_id' not in col]
    df = df.drop_duplicates(cols)
    df = create_age_group_id(df)
    df = create_sex_id(df)
    df = map_cause_cols(df)
    df = map_hosp(df)
    df.rename(columns={'drug_name': 'raw_antibiotic'}, inplace=True)


    df['resistance_old'] = df.interpretation.map({
    'Resistant': 'resistant',
    'Susceptible': 'susceptible',
    'Intermediate': 'resistant',
    'Positive': 'resistant',
    'Non Suceptible': 'resistant',
    'Negative': 'sensitive'
    }).fillna('unknown')

    df["location_id"] = 102
    df["nid"] = "NUMBER"
    df["year_id"] = df["discharge_date"].apply(lambda x: x.year)
    df['sample_id'] = (df['patient_id'].astype(str) + df['specimen_id'].astype(str))
    df['cases'] = 1

    df['hospital_type'] = 'unknown'
    df['hospital_name'] = 'unknown'

    df = get_ast_values(df)

    salmonella_list = ['salmonella_c2_serotype_newport','salmonella_i_4,[5],12:i','salmonella_cholerasuis', 'shigella_dysen_(a)',
    'salmonella_d1_serotype_javiana', 'salmonella_serotype_mississippi', 'salmonella_species_group_o','salmonella_h_nos', 
    'salmonella_kiambu', 'salmonella_i_nos', 'salmonella_species,_nos','salmonella_serotype_montevideo', 'salmonella_b_serotype_heidelberg', 'salmonella_b_serotype_st_paul', 
    'salmonella_c2_serotype_muenchen','salmonella_b_nos', 'salmonella_species_subspecies_diarizonae','salmonella_anatum', 'salmonella_c1_nos', 
    'salmonella_c2_nos', 'salmonella_serotype_oranienburg', 'salmonella_c1_serotype_infantis','salmonella_d_nos','salmonella_serotype_norwich', 
    'salmonella_c_nos', 'salmonella_e_nos','salmonella_serotype_litchfield','salmonella_serotype_havana', 'salmonella_c/d_nos', 'salmonella_serotype_java',
    'salmonella_serotype_agona', 'salmonella_b_serotype_stanley', 'salmonella_serotype_indiana', 'salmonella_uganda_serogroup_e1','salmonella_c1_serotype_braenderup', 
    'salmonella_serotype_poona', 'salmonella_c1_serotype_virchow', 'salmonella_group_iiib', 'salmonella_schwarzegrund', 'salmonella_serotype_telelkebir',
    'salmonella_e/g,_nos','salmonella_sp_serotype_agbeni','salmonella_serotype_gaminara','salmonella_b_serotype_derby','salmonella_serotype_dublin', 
    'salmonella_subspecfies_iiia','salmonella_serotype_arizonae', 'salmonella_serotype_panama','salmonella_serotype_thompson',
    'salmonella_h_type_a', 'salmonella_serotype_san_diego','salmonella_species_group_f', 'salmonella_serotype_berta', 'salmonella_inverness', 
    'salmonella_serotype_give', 'salmonella_group_z','salmonella_bareilly','salmonella_b_serotype_reading', 'salmonella_rubislaw_(serogroup_f)', 'salmonella_steftenberg', 
    'salmonella_a_nos', 'salmonella_adelaide', 'salmonella_serotype_hadar', 'salmonella_serotype_rissen','salmonella_paratyphi_c','salmonella_serotype_weltevreden', 
    'salmonella_manhattan', 'salmonella_sp_serotype_florida', 'salmonella_species_(serogroups)', 'salmonella_g_nos', 'salmonella_urbana', 
    'salmonella_serotype_monschaui', 'salmonella_group_1', 'salmonella_serotype_oslo','salmonella_bredeney,_group_b', 'salmonella_gp_2', 
    'salmonella_d1_serotype_miami', 'salmonella_serotype_ohio','salmonella_serotype_minnesota','salmonella_c1_serotype_mbandaka','salmonella_species_serotype_amunigun', 
    'salmonella_serotype_hartford', 'salmonella_serotype_london','salmonella_tallahassee', 'salmonella_serotype_tananarive',
    'salmonella_serotype_worthington', 'salmonella_lomalinda', 'salmonella_serotype_chester', 'salmonella_species_serotype_kimuenza',
    'salmonella_c2_serotype_blockley', 'salmonella_senftenberg', 'salmonella_serotype_sentfenberg', 'salmonella_species_group_r', 
    'salmonella_serotype_kentucky', 'salmonella_serogroup_nottingham','salmonella_serotype_putten','salmonella_durban', 
    'salmonella_serotype_alachua', 'salmonella_bovismorbificans', 'salmonella_e1_muenster','salmonella_serotype_fillmore','salmonella_korbol', 'salmonella_serotype_kisangani', 
    'salmonella_pensacola', 'salmonella_serotype_kampala', 'salmonella_okatie','shigella_flexneri_2a', 'salmonella_serotype_lagos', 
    'salmonella_h_type_de','salmonella_serotype_othmarschen', 'salmonella_serotype_carrau', 'salmonella_serotype_gozo','salmonella_serotype_dusseldorf',
    'salmonella_altona','salmonella_serotype_braenderup','salmonella_serotype_brandenberg','salmonella_grumpensis','salmonella_species_serotype_goldcoast','salmonella_concord','salmonella_cerro','salmonella_cotham',
    'salmonella_iv','salmonella_species_serotype_chailey','salmonella_species_group_s','salmonella_serotype_johannesburg','salmonella_hvittingfoss','salmonella_sp_serotype_uzaramo','salmonella_serotype_aberdeen',
    'salmonella_sundsvall','salmonella_clackamas','salmonella_serotype_baildon','salmonella_daytona','salmonella_serotype_soahanina',
    'salmonella_serotype_anecho','salmonella_group_3','salmonella_species_serotype_blegdam', 'salmonella_madelia', 'salmonella_group_k'
    ]

    df.raw_pathogen = df.raw_pathogen.str.lower().str.replace(' ', '_')
    df.loc[df.raw_pathogen.isin(salmonella_list), 'raw_pathogen'] = 'salmonella'

    df = create_resistance_column_clsi_2023(df)

    df['interpretation_method'] = df.apply(lambda x:'CLSI_2023' if x.resistance!='unknown' else 'other', axis=1)
    df['resistance']= df.apply(lambda x: x.resistance_old if x.resistance=='unknown' else x.resistance, axis=1)
    
    cause_columns = [x for x in list(df) if "multiple_cause" in x] + ['cause']
    demo_cols = ['nid', 'location_id', 'year_id', 'age_group_id', 'sex_id']
    biology = ['raw_specimen', 'raw_pathogen', 'raw_antibiotic', 'resistance', 'interpretation_method']
    other_values = ['hosp', 'hospital_name', 'hospital_type', 'deaths', 'sample_id', 'cases']
    lowercase_cols = ['raw_specimen', 'raw_pathogen', 'raw_antibiotic']
    resistance_testing = ['method', 'number', 'number_val', 'site_category']

    df = df[cause_columns+demo_cols + biology + other_values + resistance_testing]


    for col in lowercase_cols:
        print(col)
        df[col] = df[col].str.lower().str.strip()
    assert df.sample_id.notnull().values.all()

    df = df.drop_duplicates()



    return df

if __name__ == "__main__":
    df = clean_SOURCE()
    df['raw_pathogen']  = df.raw_pathogen.str.replace('_',' ')
    df.to_csv("FILEPATH", index=False)

    df = df.drop(['method', 'number', 'number_val', 'site_category'], axis=1)
    save_amr_data(df, phase='unmapped', ldir=L_DIR)