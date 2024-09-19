import pandas as pd
import numpy as np
import sys
import os
import argparse
from warnings import warn
import getpass
import ipdb

user = getpass.getuser()
amrpath = 'FILEPATH'.format(user)
if not amrpath in sys.path:
    sys.path.append('FILEPATH'.format(user))
from amr_prep.utils.amr_io import get_amr_data
from mcod_prep.utils.mcause_io import get_mcause_data
from datetime import datetime
from cod_prep.downloaders import locations
from cod_prep.utils import print_log_message

repo_dir = '~/FILEPATH'
MAP_DIR = repo_dir + "/FILEPATH/"
adj_dir = MAP_DIR + 'FILEPATH/'
breakpt_dir = MAP_DIR + 'FILEPATH/'

pathogen_map = pd.read_csv(f'{breakpt_dir}FILEPATH')
pathogen_map = pathogen_map.drop_duplicates('raw_pathogen')

clsi23_datasets = ['SOURCE']


def assign_original_guideline(df):

    def fit_region_to_guideline(df):


        clsi_region_ids = [100,199,9,65,5,96,159,70,134,138,124,174,120,104,167,32,192,21] 
        eucast_region_ids = [73,42,56]

        df.loc[df.region_id.isin(clsi_region_ids),'location_guideline'] = 'CLSI'
        df.loc[df.region_id.isin(eucast_region_ids),'location_guideline'] = 'EUCAST'

        return df

    df =  locations.add_location_metadata(df,['region_id','location_name'])
    df = fit_region_to_guideline(df)
    
    
    df['year_for_guidelines'] = df.year_id.apply(lambda x:x if x>=2011 else 2011 )
    df['guideline'] = df['location_guideline']+ '_'+df['year_for_guidelines'].astype(int).astype(str)

    # based on metadata
    df.loc[df.source=='SOURCE','guideline']='CLSI_XXXX'

    return df

def add_adj_factors(df):

    adj = pd.read_csv(f"{MAP_DIR}/FILEPATH")
    adj_combos = adj[['raw_antibiotic','pathogen','site_category']].drop_duplicates()

    adj_combos.rename(columns= {'pathogen':'adj_pathogen'}, inplace=True)
    
    df['site_category'] = 'other'
    df.loc[df.specimen=='urinary', 'site_category'] = 'uti'
    df.loc[df.specimen=='csf', 'site_category'] = 'csf'


    df = df.merge(pathogen_map, how = 'left', left_on='pathogen', right_on='raw_pathogen',validate='m:1')
    df = df.merge(adj_combos, how='left',left_on=['mapped_pathogen','raw_antibiotic','site_category'], right_on=['adj_pathogen','raw_antibiotic','site_category'], validate='m:1', indicator='first_match')

    # matched rows
    df_m = df[df.first_match=='both']
    # unmatched rows
    df_nm = df[df.first_match!='both']

    # preparing for additional merge, drop the column to avoid duplication 
    df_nm.drop('adj_pathogen', axis=1, inplace=True)

    # second try at merging, not based on 'pathogen_name_2'
    df_nm= df_nm.merge(adj_combos, how='left',left_on=['pathogen_name_2','raw_antibiotic','site_category'], right_on=['adj_pathogen','raw_antibiotic','site_category'], validate='m:1', indicator='second_match')
    df = pd.concat([df_m, df_nm])

    df.drop(['raw_pathogen','mapped_pathogen','pathogen_name_2', 'first_match','second_match'],axis=1, inplace=True)

    df['interpretation_method'] = df.interpretation_method.str.replace(' ','_').str.upper()

    df = assign_original_guideline(df)

    # assign the intepretation_method as the guideline if it is a speicific one, not "other"
    df.loc[(df.interpretation_method!='other')&(~df.interpretation_method.isna())&(df.interpretation_method!='OTHER'), 'guideline'] = df.interpretation_method

    
    # keep the guideline names consistant 
    df['guideline'] = df.guideline.str.upper()
    adj['guideline'] = adj.guideline.str.replace(' ', '_')

    adj = adj.drop_duplicates()
    df = df.merge(adj,left_on=['adj_pathogen','raw_antibiotic','site_category','guideline'], right_on=['pathogen','raw_antibiotic','site_category','guideline'], how='left', validate='m:1')

    df.loc[df.adj_factor.isna(), 'adj_factor'] = 1 


    df.drop(['pathogen_y','datapoints','method','adj_pathogen'],axis=1,inplace=True)
    df.rename(columns={'pathogen_x':'pathogen'},inplace=True)

    return df

def testing(df, n):
    
    assert len(df.guideline.unique()) == 25 
    assert len (df.source.unique()) == n


def check_adj_factors(check):
    '''
    check that that adj factors make sense
    '''
    adj = pd.read_csv(f"{MAP_DIR}/FILEPATH")
    g = g.replace(' ',' _')

    # check that all adj factors used are factors present in the relevant guideline 
    for g in adj.guideline.unique():
        g = g.replace(' ',' _')
        for a in check.loc[check.guideline==g,'adj_factor'].unique():
            assert a in adj.loc[adj.guideline==g, 'adj_factor'].unique()

    print_log_message('checking adjustments')  
    for i in range(len(adj)):
    
        p = adj.iloc[i, 0]
        ab = adj.iloc[i, 1]
        sc = adj.iloc[i, 2]
        g = adj.iloc[i, 3]

        adj_f = adj.iloc[i, 4]

        print(f"{p}, {ab}, {sc}, {g}")
        try:
            assert check.loc[(check.pathogen==p)&(check.raw_antibiotic==ab)&(check.site_category==sc)&(check.guideline==g),'adj_factor'].unique()[0] == adj_f
               
        except:
            if adj_f!=1:
                psb_p = check.loc[(check.raw_antibiotic==ab)&
                                  (check.site_category==sc)&(check.guideline==g)&(check.adj_factor==adj_f), 'pathogen'].unique()
                print(f'this pathogen name not present in the data, these were adjuted with the same adj factor {psb_p}')

            continue
        

        for g in adj.guideline.unique():
            if g!='CLSI_XXXX':
                    adjs = check.loc[check.guideline==g, 'adj_factor'].unique()
                    assert len(adjs)>1, f'only one adj facotor {adj[0]} is available for the guideline {g}'




def calculate_ressum(resdat, make_guideline_adj = True):

    print_log_message('calculating sums of resistance...')
    
    resdat.loc[resdat.interpretation_method.isna(), 'interpretation_method'] = 'other'
    
    if make_guideline_adj:
   
        ressum_adj = resdat.groupby(['location_id', 'year_id','nid','source','specimen', 'pathogen', 'raw_antibiotic','abx_class','interpretation_method','hospital_type'], group_keys=False, as_index = False).apply(

            lambda x: pd.Series({'resistant':sum(x.loc[x['resistance'] == 'resistant', 'cases']),
                                'susceptible':sum(x.loc[x['resistance'] == 'susceptible', 'cases']),
                                'cases':sum(x['cases'])})).reset_index()

        ressum_adj = add_adj_factors(ressum_adj)    

        ressum_adj['resistance_fraction'] = ressum_adj['resistant']/ressum_adj['cases']
        ressum_adj['resistance_fraction_adj'] = ressum_adj['resistance_fraction']*ressum_adj['adj_factor']

        ressum_adj.loc[ressum_adj.resistance_fraction_adj>1, 'resistance_fraction_adj'] = 0.99

        check_adj_factors(ressum_adj)
        testing(ressum_adj)


        ressum_adj['resistant_adj'] = ressum_adj['cases']*ressum_adj['resistance_fraction_adj']
        ressum_adj['susceptible_adj'] = ressum_adj['cases'] - ressum_adj['resistant_adj']

        ressum_adj['diff'] = ressum_adj['resistant']-ressum_adj['resistant_adj']

        ressum_adj.drop('raw_antibiotic',axis=1, inplace=True)

        keep_cols = [ 'location_id', 'year_id', 'nid','source', 'specimen', 'pathogen',
               'abx_class','hospital_type',  'cases',
                'resistant_adj',
               'susceptible_adj']

        ressum_adj = ressum_adj[keep_cols].groupby(['location_id', 'year_id', 'nid', 'source', 'pathogen',
               'abx_class','hospital_type']).sum().reset_index()

        ressum_adj.rename(columns={'resistant_adj':'resistant', 'susceptible_adj':'susceptible'}, inplace=True)
        
        output_df =  ressum_adj
    
    else:
        
        ressum = resdat.groupby(['location_id', 'year_id', 'nid', 'source', 'pathogen', 'abx_class','hospital_type'], group_keys=False, as_index = False).apply(
        lambda x: pd.Series({'resistant':sum(x.loc[x['resistance'] == 'resistant', 'cases']),
                        'susceptible':sum(x.loc[x['resistance'] == 'susceptible', 'cases']),
                        'cases':sum(x['cases'])})).reset_index()
        output_df = ressum

    return output_df

def set_interpretation_method(df, clsi23_datasets):

    '''create and update interpretation_method column in datasets that are missing it'''
    if 'interpretation_method' in df.columns:
        df.loc[df.interpretation_method.isna(), 'interpretation_method'] = 'other'
        
    if 'interpretation_method' not in df.columns:
            
        df.loc[df.source.isin(clsi23_datasets),'interpretation_method']='CLSI_2023'
        df.loc[~df.source.isin(clsi23_datasets),'interpretation_method']='other'
            
        print(f"{df.loc[df.interpretation_method.isna(), 'source'].unique()}  were interperted with CLSI_2023" )
        df.loc[(df.interpretation_method.isna())&(df.source.isin(clsi23_datasets)),'interpretation_method']='CLSI_2023'
        
        print(f"{df.loc[df.interpretation_method.isna(),'source'].unique()} were assigned interpretation" )
        df.loc[df.interpretation_method.isna(),'interpretation_method']='other'
        
    return df

print_log_message('Pulling amr data...')
amr_dat = get_amr_data(redistribute_pathogens=False)
amr_dat = set_interpretation_method(amr_dat, clsi23_datasets)

n_datasets=len(amr_dat.source.unique())

for s in ['citrobacter','serratia','providencia','proteus','enterobacter','moraxella','morganella','shigella']:
    amr_dat.loc[amr_dat.pathogen.str.contains(s), 'pathogen'] = s + '_spp'

for s in ['campylobacter','listeria']:
    amr_dat.loc[amr_dat.pathogen.str.contains(s), 'pathogen'] = s

for s in ['group_b','group_a']:
    amr_dat.loc[amr_dat.pathogen.str.contains(s), 'pathogen'] = s + '_strep'
    
amr_dat.loc[amr_dat.pathogen.str.contains('enterococcus') & 
            ~amr_dat.pathogen.isin(['enterococcus_faecium','enterococcus_faecalis']), 'pathogen'] = 'enterococcus_spp'

amr_dat.loc[amr_dat.pathogen.str.contains('gonorrheae'), 'pathogen'] = 'neisseria_gonorrhoeae'
amr_dat.loc[amr_dat.pathogen.str.contains('baumanii'), 'pathogen'] = 'acinetobacter_baumannii'


# overlap with other sources
amr_dat = amr_dat.loc[~amr_dat.source.isin(['SOURCE']),]
amr_dat.loc[(amr_dat['pathogen'] == 'acinetobacter_spp') & (amr_dat['source'] == 'SOURCE'), 'pathogen'] = 'acinetobacter_baumanii'

amr_dat = amr_dat.loc[~((amr_dat.source=='SOURCE') & \
                        (amr_dat.pathogen == 'staphylococcus_aureus') & \
                        (amr_dat.abx_class == 'vancomycin')),]

resdat = amr_dat.loc[(amr_dat['resistance'].notna()) & (amr_dat['resistance'].isin(['sensitive', 'resistant', 'susceptible'])), :]

del amr_dat

resdat['resistance'] = resdat['resistance'].str.replace('sensitive', 'susceptible')

if resdat['hospital_type'].isna().values.any():
    warn(f"{resdat.loc[resdat['hospital_type'].isna(), 'source']} have unknown hospital type")
if resdat['specimen'].isna().values.any():
    warn(f"{resdat.loc[resdat['specimen'].isna(), 'source']} have unknown specimen type")

resdat.loc[(resdat['hospital_type'].isna() & resdat.source.isin(['Sentry','pfizer'])),'hospital_type'] = 'mixed'
resdat.loc[resdat['hospital_type'].isna(),'hospital_type'] = 'unknown'

resdat.loc[resdat['specimen'].isna(),'specimen'] = 'unknown'


def remove_intrinsic_resistance(df):
    print_log_message('Removing intrinsic resistance')
    
    ir  = pd.read_csv(f"{MAP_DIR}FILEPATH")
    ir_drug = ir[ir.type=='drug']
    ir_class = ir[ir.type=='class']
    
    init_len = len(df)
    
    cols = df.columns
    df = df.merge(ir_drug, how='left', left_on=['pathogen','raw_antibiotic'], right_on=['pathogen','antibiotic'], validate='m:1', indicator=True)
    
    df = df[df._merge!='both']
    
    # return to the same cols as before merge

    df = df[cols]
    
    df = df.merge(ir_class, how='left', left_on=['pathogen','abx_class'], right_on=['pathogen','antibiotic'], validate='m:1', indicator=True)
    df = df[df._merge!='both']
    
    df = df[cols]
    
    print_log_message(f"{len(df)-init_len} rows were removed")
    
    return df

resdat = remove_intrinsic_resistance(resdat)
ressum = calculate_ressum(resdat=resdat, make_guideline_adj=True)

needrespath = ['acinetobacter_baumanii',
		   'acinetobacter_baumannii',
              'campylobacter',
              'chlamydia_spp',
              'citrobacter_spp',
              'clostridium_difficile',
              'cryptosporidiosis',
              'enterobacter_spp',
              'enterococcus_faecium',
              'enterococcus_faecalis',
              'enterococcus_spp',
              'escherichia_coli',
              'group_a_strep',
              'group_b_strep',
              'haemophilus_influenzae',
              'klebsiella_pneumoniae',
              'listeria',
              'moraxella_spp',
              'mycoplasma',
              'neisseria_meningitidis',
              'morganella_spp',
              'neisseria_gonorrheae','neisseria_gonorrhoeae',
              'proteus_spp',
              'providencia_spp',
              'pseudomonas_aeruginosa',
              'pseudomonas_spp',
              'serratia_spp',
              'legionella_spp',
              'non_typhoidal_salmonellae',
              'salmonella_typhi',
              'salmonella_paratyphi',
              'shigella_spp',
              'streptococcus_pneumoniae',
              'staphylococcus_aureus',
              'mycobacterium_tuberculosis',
              'salmonella_ints_non_typhi/paratyphi']

ressum2 = ressum.loc[ressum['pathogen'].isin(needrespath),:].reset_index()

ressum_old = pd.read_csv('FILEPATH')
ressum_old.to_csv(f'FILEPATH', index=False)

ressum2.to_csv('FILEPATH', index = False)

