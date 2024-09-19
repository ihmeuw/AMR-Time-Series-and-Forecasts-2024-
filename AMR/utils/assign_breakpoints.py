
import pandas as pd
import numpy as np
from cod_prep.utils import print_log_message
from pathlib import Path


current_fp =  "FILEPATH"
root_dir =  "FILEPATH"
repo_dir =  "FILEPATH"
MAP_DIR =  "FILEPATH"
L_DIR =  "FILEPATH"
breakpt_dir =  "FILEPATH"
intermediate_dir = "FILEPATH"
raw_dir =  "FILEPATH"

breakpt_map = pd.read_csv("FILEPATH")
pathogen_map = pd.read_csv("FILEPATH")
pathogen_map = pathogen_map.drop_duplicates('raw_pathogen')
ab_brkpt_list = breakpt_map.antibiotic.unique()


def print_out_antibiotic_names(ab_brkpt = ab_brkpt_list):

    ab_brkpt_list = list(ab_brkpt)
    ab_brkpt_list.sort()

    print(ab_brkpt_list)

def validate_format_antibiotic_names(df, ab_brkpt = ab_brkpt_list):
    ab_names_dict = {'co-amoxiclav':'amoxicillin/clavulanic acid',
                    'amoxicillin-clavulanic acid':'amoxicillin/clavulanic acid',
                    'amoxicillin clavulanate':'amoxicillin/clavulanic acid',
                    'amoxicillin /clavulanic acid':'amoxicillin/clavulanic acid',
                    'amoxicillin/ca':'amoxicillin/clavulanic acid',
                    'high level streptomycin':'streptomycin-high',
                    'streptomycin/high':'streptomycin-high',   
                    'cefuroxime axetil':'cefuroxime',
                    'penicillin g': 'benzylpenicillin',
                    'piptaz':'piperacillin/tazobactam',
                    'tazocin':'piperacillin/tazobactam',
                    'daptomycin e ttest':'daptomycin',
                    'ampicillin sulbactam':'ampicillin/sulbactam',
                    'piperacillin tazobactam': 'piperacillin/tazobactam',
                    'gentamicin/high':'gentamicin-high',
                    'cephalotin':'cephalothin',
                    'augmentin':'amoxicillin/clavulanic acid',
                'ceftazidime clsi':'ceftazidime',
                'nalidixic':'nalidixic acid',
               'trimethorpim' :'trimethoprim',
                'sulfonamide':'sulfonamide',
                'gentamycin':'gentamicin',
                'netilmicin clsi':'netilmicin',
                'cefotaxim clsi':'cefotaxime',
                'sulphonamides.1':'sulfonamide',
                'suphamethoxazole':'sulfamethoxazole',
                'azithromycine':'azithromycin',
                'sulphonamides':'sulfonamide',
                'high level gentamici':'gentamicin-high',
                'imp carbapenemase':'imipenem',
                'amoxycillin clavulanate':'amoxicillin/clavulanic acid',
                'meropenem vaborbactam':'meropenem/vaborbactam',
                'ceftazidime avibactam':'ceftazidime/avibactam',
                'quinupristin dalfopristin':'quinupristin/dalfopristin',
                'trimethoprim sulfa':'trimethoprim/sulfamethoxazole',
                'cotrimoxazole':'trimethoprim/sulfamethoxazole',
                'trimethoprim/sulfa (sxt)':'trimethoprim/sulfamethoxazole',
                'amp/sulbactam':'ampicillin/sulbactam',
                'streptomycin high/level':'streptomycin-high',
                'streptomycin 1000':'streptomycin-high'
                }

    
    df['raw_antibiotic'] = df['raw_antibiotic'].str.strip().str.lower()
    df['raw_antibiotic'] = df['raw_antibiotic'].apply(lambda x: ab_names_dict[x] if x in ab_names_dict.keys() else x)
    df['raw_antibiotic'] = df['raw_antibiotic'].str.replace('-','/')

    ab_names_dict2 = {'streptomycin/high': 'streptomycin-high',
                    'gentamicin/high':'gentamicin-high', 
                    'streptomycin high/level':'streptomycin-high'}    
    df['raw_antibiotic'] = df['raw_antibiotic'].apply(lambda x: ab_names_dict2[x] if x in ab_names_dict2.keys() else x)

    unmatched_ab = []
    for i in df.raw_antibiotic.unique():
        if i not in ab_brkpt:
            unmatched_ab.append(i)
    if len(unmatched_ab)>0:
        print (f'The following antibiotic names were not matched with breakpoint antibiotics: {unmatched_ab}')

    return df

def missing_pathogen_names(df, pathogen_map=pathogen_map):
    not_present = []
    for i in df.raw_pathogen.unique():
        if i not in pathogen_map.raw_pathogen.unique():
            not_present.append(i)
    print(f'The following pathogen names were not matched:\n{not_present}')



def assign_breakpoints(df, g, breakpt_map=breakpt_map, pathogen_map=pathogen_map):

    if 'disk_dose' not in df.columns:
        df['disk_dose'] = 'not_disk'
    df.loc[(df.disk_dose.isna())&(df.method=='MIC'), 'disk_dose']='not_disk'

    breakpt_map = breakpt_map[breakpt_map.use_breakpoint==1]
    assert g in breakpt_map.guideline.unique(), f'{g} is not in the guideline list. Please provide one of the following guidelines: {breakpt_map.guideline.unique()}'
    b_map = breakpt_map[['guideline', 'method','site_category',
           'antibiotic', 'breakpoint_S', 
           'pathogen', 'use_breakpoint', 'disk_dose_1', 'disk_dose_2']]

    b_map = b_map[b_map.guideline==g]
    b_map.loc[b_map.method=='MIC', 'disk_dose_1'] = 'not_disk'

    print_log_message('formatting pathogen names')
    df['raw_pathogen'] = df['raw_pathogen'].str.strip().str.lower()
    df['raw_pathogen'] = df['raw_pathogen'].str.replace(' ','_').str.lower()

    missing_pathogen_names(df)

    df = df.merge(pathogen_map, how='left', validate="m:1")

    df_a = df.merge(b_map, how='left', left_on = ['method','site_category','mapped_pathogen','raw_antibiotic','disk_dose'],
             right_on=[ 'method','site_category','pathogen','antibiotic','disk_dose_1'],validate='m:1',
            indicator=True)

    matched_a = df_a[df_a._merge=='both']
   
    not_matched = df_a[df_a._merge=='left_only']

    cols_to_remove = ['antibiotic', 'guideline','breakpoint_S','pathogen', 'use_breakpoint', 'disk_dose_1', 'disk_dose_2','_merge']
    not_matched.drop(cols_to_remove, axis=1, inplace=True)

    df_b = not_matched.merge(b_map, how='left', left_on = ['method','site_category','pathogen_name_2','raw_antibiotic','disk_dose'],
         right_on=[ 'method','site_category','pathogen','antibiotic','disk_dose_1'],validate='m:1', indicator=True)

    matched_b = df_b[df_b._merge=='both']
    not_matched = df_b[df_b._merge=='left_only']

    not_matched['site_category'] = 'other'
    not_matched.drop(cols_to_remove, axis=1, inplace=True)

    df_c = not_matched.merge(b_map, how='left', left_on = ['method','site_category','mapped_pathogen','raw_antibiotic','disk_dose'],
         right_on=[ 'method','site_category','pathogen','antibiotic','disk_dose_1'], validate='m:1',
        indicator=True)

    matched_c = df_c[df_c._merge == 'both']
    not_matched = df_c[df_c._merge == 'left_only']

    not_matched.drop(cols_to_remove, axis=1, inplace=True)

    df_d = not_matched.merge(b_map, how='left', left_on = ['method','site_category','pathogen_name_2','raw_antibiotic','disk_dose'],
         right_on=[ 'method','site_category','pathogen','antibiotic','disk_dose_1'],
        indicator=True, validate='m:1')

    matched_d = df_d[df_d._merge=='both']
    not_matched = df_d[df_d._merge=='left_only']

    len(df) == (len(matched_a)+len(matched_b)+len(matched_c)+len(matched_d)+len(not_matched))
    df = pd.concat([matched_a, matched_b, matched_c, matched_d, not_matched])
    percent_matched = len(df[~df.breakpoint_S.isna()])/len(df)*100
    print_log_message(f'{percent_matched}% of dataset were assigned breakpoints')

    return df











