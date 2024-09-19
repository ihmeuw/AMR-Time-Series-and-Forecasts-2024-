import glob
from pathlib import Path

import importlib
import db_queries.api.public as db

import re
from cod_prep.utils.formatting import ages


import os
import pandas as pd
import numpy as np
from cod_prep.claude.configurator import Configurator
from cod_prep.utils import report_if_merge_fail, print_log_message, cod_timestamp

from amr_prep.utils.amr_io import save_amr_data

repo_dir = "FILEPATH"
MAP_DIR =  "FILEPATH"
L_DIR =  "FILEPATH"
breakpt_dir =  "FILEPATH"
intermediate_dir =  "FILEPATH"
raw_dir =  "FILEPATH"


def  create_breakpoint_map():
    
    breakpt_map = pd.read_excel( "FILEPATH")

    breakpt_map['ab_name'] = breakpt_map.ab_name.str.strip().str.lower() 
    
    breakpt_map['pathogen'] = breakpt_map.mo_name.str.strip().str.lower() 
    breakpt_map['pathogen'] = breakpt_map.pathogen.str.replace(' ','_')

    breakpt_map['site_category'] = 'other'

    breakpt_map.loc[breakpt_map.uti==True, 'site_category']='uti'
    breakpt_map.loc[breakpt_map.site=='Meningitis', 'site_category'] = 'csf' 

    breakpt_map = adjust_brkpt_map_v2(breakpt_map)

    breakpt_map = breakpt_map[breakpt_map.type=='human']

    breakpt_map[['disk_dose_1', 'disk_dose_2']] = breakpt_map.disk_dose.str.split('/', expand=True)
    for col in ['disk_dose_1', 'disk_dose_2']:
        breakpt_map[col] = breakpt_map[col].str.replace('ug','')
        
        breakpt_map[col] = breakpt_map[col].str.replace(' units','')
        breakpt_map[col] = breakpt_map[col].str.replace(' unit','')
        breakpt_map[col] = breakpt_map[col].str.replace('10-15','10')
        breakpt_map[col] = breakpt_map[col].astype(float)

    breakpt_map.rename(columns={'ab_name':'antibiotic'}, inplace=True)
    keep_cols = ['guideline', 'type', 'method', 'site', 'site_category','mo_name', 
           'ab', 'antibiotic','disk_dose', 'breakpoint_S', 'breakpoint_R',
           'uti','pathogen', 'use_breakpoint', 'disk_dose_1', 'disk_dose_2']


    breakpt_map = breakpt_map[keep_cols]
    
    print_log_message('Saving breakpoint map')
    breakpt_map.to_csv( "FILEPATH", index=False)


def adjust_brkpt_map_v2(breakpt_map):
    
    gdls = ['CLSI 2011', 'CLSI 2012', 'CLSI 2013',
       'CLSI 2014', 'CLSI 2015', 'CLSI 2016', 'CLSI 2017', 'CLSI 2018',
       'CLSI 2019', 'CLSI 2020', 'CLSI 2021', 'CLSI 2022', 'CLSI 2023',
       'EUCAST 2012', 'EUCAST 2013', 'EUCAST 2014', 'EUCAST 2015',
       'EUCAST 2016', 'EUCAST 2017', 'EUCAST 2018', 'EUCAST 2019',
       'EUCAST 2020', 'EUCAST 2021', 'EUCAST 2022']


    clsi_gdls = [g for g in gdls if re.match(r'^CLSI', g)]
    eucast_gdls = [g for g in gdls if re.match(r'^EUCAST', g)]
    
    breakpt_map['use_breakpoint'] = 1

    breakpt_map.loc[(breakpt_map.guideline.isin(['CLSI 2023']))
                &(breakpt_map.pathogen=='mycobacterium_avium-intracellulare_complex')
                &(breakpt_map.ab_name=='amikacin')
               &(breakpt_map.site=='Liposomal, Inhaled'), 'use_breakpoint'] = 0

    
    breakpt_map.loc[(breakpt_map.guideline.isin(clsi_gdls))
                &(breakpt_map.pathogen=='streptococcus_pneumoniae')
                &(breakpt_map.ab_name=='benzylpenicillin')
               &(breakpt_map.site=='Oral'), 'use_breakpoint'] = 0

    
    breakpt_map.loc[(breakpt_map.guideline.isin(clsi_gdls))
                &(breakpt_map.pathogen=='streptococcus_pneumoniae')
                &(breakpt_map.ab_name=='cefuroxime')
               &(breakpt_map.site=='Oral'), 'use_breakpoint'] = 0
    
    breakpt_map.loc[(breakpt_map.guideline.isin(['CLSI 2020', 'CLSI 2021', 'CLSI 2022','CLSI 2023']))
                &(breakpt_map.pathogen=='enterobacterales')
                &(breakpt_map.ab_name=='cefuroxime')
               &(breakpt_map.site=='Oral'), 'use_breakpoint'] = 0

    breakpt_map.loc[(breakpt_map.guideline=='CLSI 2012')
                &(breakpt_map.pathogen=='staphylococcus')
                &(breakpt_map.ab_name=='cefuroxime')
               &(breakpt_map.site=='Intravenous'),'use_breakpoint'] = 0
    
    breakpt_map.loc[(breakpt_map.guideline.isin(['CLSI 2012','CLSI 2013']))
                &(breakpt_map.pathogen=='salmonella')
                &(breakpt_map.ab_name=='ciprofloxacin')
               &(breakpt_map.site.isin(['Intestinal','Extraintestinal'])), 'use_breakpoint']=0

    breakpt_map.loc[(breakpt_map.guideline.isin(eucast_gdls))&(breakpt_map.pathogen=='haemophilus')
                &(breakpt_map.ab_name=='delafloxacin')
               &(breakpt_map.site!='Pneumonia'),'use_breakpoint'] = 0

    breakpt_map.loc[(breakpt_map.guideline.isin(eucast_gdls))&(breakpt_map.pathogen=='haemophilus')
                    &(breakpt_map.ab_name=='amoxicillin/clavulanic acid')
                   &(breakpt_map.site=='IV'),'use_breakpoint'] = 0

    breakpt_map.loc[(breakpt_map.guideline.isin(eucast_gdls))&(breakpt_map.pathogen=='haemophilus')
                    &(breakpt_map.ab_name=='amoxicillin')
                   &(breakpt_map.site.isin(['IV','Non-meningitis'])),'use_breakpoint'] = 0

    breakpt_map.loc[(breakpt_map.guideline.isin(eucast_gdls))&(breakpt_map.pathogen=='haemophilus')
                    &(breakpt_map.ab_name=='cefuroxime')
                   &(breakpt_map.site=='IV'),'use_breakpoint'] = 0
                    
    breakpt_map.loc[(breakpt_map.guideline.isin(eucast_gdls))&(breakpt_map.pathogen=='moraxella_catarrhalis')
                    &(breakpt_map.ab_name=='cefuroxime')
                   &(breakpt_map.site=='IV'),'use_breakpoint'] = 0
    

    breakpt_map.loc[(breakpt_map.guideline.isin(eucast_gdls))&(breakpt_map.pathogen=='haemophilus_influenzae')
                    &(breakpt_map.ab_name=='cefuroxime')
                   &(breakpt_map.site=='IV'),'use_breakpoint'] = 0
    
    breakpt_map.loc[(breakpt_map.guideline.isin(eucast_gdls))&(breakpt_map.pathogen=='streptococcus_pneumoniae')
                    &(breakpt_map.ab_name=='cefuroxime')
                   &(breakpt_map.site=='IV'),'use_breakpoint'] = 0


    breakpt_map.loc[(breakpt_map.guideline=='EUCAST 2013')&(breakpt_map.pathogen=='streptococcus_pneumoniae')
                    &(breakpt_map.ab_name=='cefuroxime')
                   &(breakpt_map.site.isna()),'use_breakpoint'] = 0


    breakpt_map.loc[(breakpt_map.guideline=='EUCAST 2023')
                &(breakpt_map.ab_name=='amoxicillin/clavulanic acid')
                &(breakpt_map.pathogen=='enterobacterales')
               &(breakpt_map.method=='MIC')
                &(breakpt_map.site=='Intravenous'), 'use_breakpoint'] = 0


    return breakpt_map