
# The goal of this script is to filter the intrinsic resistance combinations in the CLSI guideline for the one we use with our combinations
# and adjutst the name
# Anna 

import re
import os
import glob

import pandas as pd
import numpy as np
from cod_prep.utils import print_log_message
import matplotlib.pyplot as plt
import seaborn as sns
from cod_prep.downloaders import locations

repo_dir = '~/amr'
MAP_DIR = repo_dir + "/maps/"
L_DIR = "/ihme/limited_use/LIMITED_USE/LU_AMR/"
breakpt_dir = MAP_DIR + 'resistance_cutoff_mapping/'
intermediate_dir = breakpt_dir + 'intermediate_data/'
raw_dir = breakpt_dir + 'raw_data/'

breakpt_map = pd.read_csv(f'{breakpt_dir}breakpoint_map_2023.csv')
p_map_b = pd.read_csv(f'{breakpt_dir}pathogen_name_map_for_breakpoints.csv')

output_dir = '~/notebook_outputs/'

ab_map = pd.read_csv(f'{MAP_DIR}antibiotic_map.csv')
p_map = pd.read_csv(f"{MAP_DIR}pathogen_map.csv")

# filter for only the combinations we use
combos = pd.read_csv(f"{MAP_DIR}bug_drug_combos.csv", encoding='latin-1')
combos = combos[combos.kevin_catrin_exclude==0]
combos  = combos[['pathogen','abx_class']]

p_map['raw_pathogen'] = p_map.raw_pathogen.str.replace(' ','_')

int_resis = pd.read_excel(f"{breakpt_dir}intrinsic_resistant.xlsx")
int_resis['pathogen'] = int_resis.microorganism.str.lower().str.replace(' ', '_')
int_resis['raw_antibiotic'] = int_resis.antibiotic.str.lower()
int_resis['mapped_pathogen'] = int_resis['pathogen']

int_resis.loc[int_resis.pathogen=='proteus','mapped_pathogen']  = 'proteus_spp'
int_resis.loc[int_resis.pathogen=='shigella','mapped_pathogen'] = 'shigella_spp'
int_resis.loc[int_resis.pathogen=='serratia','mapped_pathogen'] = 'serratia_spp'
int_resis.loc[int_resis.pathogen=='citrobacter', 'mapped_pathogen'] = 'citrobacter_spp'
int_resis.loc[int_resis.pathogen=='enterobacter', 'mapped_pathogen'] = 'enterobacter_spp'
int_resis.loc[int_resis.pathogen=='enterococcus', 'mapped_pathogen'] = 'enterococcus_spp'

order_list = ['citrobacter','enterobacter','enterococcus','serratia','proteus','shigella']

for o in order_list:
    print(o)
    for i in int_resis.loc[int_resis.pathogen.str.contains(o), 'pathogen'].unique():
        if i in combos.pathogen.unique():
            print(i)
            continue
        int_resis.loc[int_resis.pathogen==i, 'mapped_pathogen'] = f"{o}_spp"

# group a: 'streptococcus_pyogenes'
# group b: 'streptococcus_agalactiae'

int_resis.loc[int_resis.pathogen.str.contains('agalactiae'), 'pathogen'].unique()
int_resis.loc[int_resis.pathogen.str.contains('streptococcus_group_a'), 'mapped_pathogen'] = 'group_a_strep'
int_resis.loc[int_resis.pathogen.str.contains('streptococcus_group_b'), 'mapped_pathogen'] = 'group_b_strep'
int_resis.loc[int_resis.pathogen.str.contains('streptococcus_pyogenes'), 'mapped_pathogen'] = 'group_a_strep'
int_resis.loc[int_resis.pathogen.str.contains('streptococcus_agalactiae'), 'mapped_pathogen'] = 'group_b_strep'

nt_sal = p_map.loc[p_map.pathogen=='non_typhoidal_salmonellae','raw_pathogen'].unique()
pt_sal = p_map.loc[p_map.pathogen=='salmonella_paratyphi','raw_pathogen'].unique()
t_sal = p_map.loc[p_map.pathogen=='salmonella_typhi','raw_pathogen'].unique()



int_resis.loc[int_resis.pathogen.isin(nt_sal)]

int_resis.loc[int_resis.pathogen.isin(nt_sal), 'mapped_pathogen'] = 'non_typhoidal_salmonellae'
int_resis.loc[int_resis.pathogen.isin(pt_sal), 'mapped_pathogen'] = 'salmonella_paratyphi'
int_resis.loc[int_resis.pathogen.isin(t_sal), 'mapped_pathogen'] = 'salmonella_typhi'

int_resis = int_resis.merge(ab_map, how='left')

int_resis = int_resis.merge(pathogen_map, how='left',left_on='pathogen', right_on='raw_pathogen')

int_resis_1 = int_resis.merge(combos, left_on=['mapped_pathogen','abx_class'], right_on=['pathogen','abx_class'], how='inner')

int_resis_1 = int_resis_1[['pathogen_x','raw_antibiotic','mapped_pathogen','abx_class']].drop_duplicates().rename(columns={'pathogen_x':'pathogen'})

int_resis_1.to_csv(f"{output_dir}intrinsic_resistance_filtered.csv", index=False)