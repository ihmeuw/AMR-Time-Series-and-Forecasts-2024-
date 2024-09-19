import re
import pandas as pd
import numpy as np
from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import get_current_location_hierarchy
from cod_prep.utils import print_log_message
from cod_prep.utils.formatting import ages
from amr_prep.utils.amr_io import save_amr_data
from amr_prep.utils import assign_breakpoints

CONF = Configurator()
SOURCE = "SOURCE"
L_DIR = "FILEPATH"
repo_dir = "FILEPATH"
MAP_DIR = repo_dir + "FILEPATH"

breakpt_dir = "FILEPATH"
intermediate_dir = "FILEPATH"
raw_dir = "FILEPATH"


def read_in_data():
    clinic = pd.read_csv("FILEPATH")
    gnb = pd.read_excel("FILEPATH")
    gpb = pd.read_excel("FILEPATH")

    clinic = clinic.drop_duplicates()
    gnb = gnb.drop_duplicates()
    gpb = gpb.drop_duplicates()

    allg = pd.concat([gnb, gpb])

    clinic.rename(columns={'Anon_BabyID': 'Anon_baby'}, inplace=True)
    for df in [clinic, allg]:
        df['Anon_baby'] = df['Anon_baby'].astype(str)

    df = allg.merge(clinic, how='left', on='Anon_baby')

    return df


def format_age_group_id(df):
    df['age'] = (pd.to_datetime(df['clindate1']) - pd.to_datetime(df['babydob'])).dt.days
    df.loc[df['age'].isna(), 'age'] = df['ageatoutcomedays']

    df['age_unit'] = 'day'

    age_formatter = ages.PointAgeFormatter()
    df = age_formatter.run(df)

    df['age_group_id'] = df['age_group_id'].astype(int)
    df.loc[df['age'].isna(), 'age_group_id'] = 390

    print_log_message(
        f"{(len(df.loc[df['age_group_id'] == 390]) / len(df)) * 100}% of rows are missing age")

    return df


def format_sex_id(df):
    df.rename(columns={'gender': 'sex_id'}, inplace=True)

    df.loc[df['sex_id'].isna(), 'sex_id'] = 9

    print_log_message(
        f"{(len(df.loc[df['sex_id'] == 9]) / len(df)) * 100}% of rows are missing sex")

    return df


def format_location_id(df):
    site_to_gbd_loc_dict = {
        'BC': 'Bangladesh',
        'BK': 'Bangladesh',
        'ES': 'Addis Ababa',
        'IN': 'West Bengal',
        'NK': 'Kano',
        'NN': 'FCT (Abuja)',
        'NW': 'FCT (Abuja)',
        'PC': 'Islamabad Capital Territory',
        'PP': 'Islamabad Capital Territory',
        'RK': 'Rwanda',
        'RU': 'Rwanda',
        'ZAT': 'Western Cape',
    }

    df['location_name'] = df['Site'].map(site_to_gbd_loc_dict)

    lh = get_current_location_hierarchy()

    df = df.merge(lh[['location_id', 'location_name']], how='left', on='location_name')

    assert df['location_id'].notnull().values.all()

    return df


def mark_deaths(df):
    df.loc[df['outcome'] == 1, 'deaths'] = 0
    df.loc[df['outcome'] == 2, 'deaths'] = 1
    df.loc[df['outcome'].isin([3, np.nan]), 'deaths'] = np.nan

    return df


def create_cause_cols(df):

    df.loc[df['clinicassum'] == 1, 'cause'] = 'sepsis'
    df.loc[df['zp'] == 1, 'multiple_cause_1'] = 'perinatal asphyxia'
    df.loc[df['pprom'] == 2, 'multiple_cause_2'] = 'preterm'

    for col in ['cause', 'multiple_cause_1', 'multiple_cause_2']:
        df.loc[df[col].isna(), col] = 'none'

    return df


def create_specimen_pathogen_col(df):
    df.rename(columns={'FinalSpeciesID': 'raw_pathogen'}, inplace=True)

    df['raw_specimen'] = 'blood culture'

    return df



def create_hosp(df):
    df['day_diff'] = (pd.to_datetime(df['clindate1']) -
                      pd.to_datetime(df['dateofadmissionintohospit'])).dt.days

    df.loc[df['day_diff'] > 2, 'hosp'] = 'community'
    df.loc[df['day_diff'] <= 2, 'hosp'] = 'hospital'
    df.loc[df['hosp'].isna(), 'hosp'] = 'unknown'

    return df


def create_nid(df):
    nids = pd.read_excel("FILEPATH")
    nids['nid'] = nids['NID']
    nids['title'] = nids.Title.str.split(' ')
    nids['site'] = nids.title.str[2]
    nids.loc[nids.site == 'Burden', 'site'] = nids['title'].str[1]
    nids['year_id'] = nids.title.str[-1].astype(int)
    site_to_full_name = {
        'BC': 'Chittagong',
        'BK': 'Kumudini',
        'ES': 'Addis',
        'IN': 'Kolkata',
        'NK': 'Kano',
        'NN': 'Abuja',
        'NW': 'Abuja',
        'PC': 'Islamabad',
        'PP': 'Islamabad',
        'RK': 'Muhanga',
        'RU': 'Kigali',
        'ZAT': 'Cape',
    }
    df['site'] = df['Site'].map(site_to_full_name)
    assert df['site'].notnull().all()
    df = df.merge(
        nids[['site', 'year_id', 'nid']],
        how='left', on=['year_id', 'site'], validate='many_to_one'
    )
    assert df.nid.notnull().all()
    return df


def add_hospital_type(df): 
    hosp_dictionary = { 
        'BC': 'Chittagong Ma o Shishu Hospital', 
        'BK': "Kumudini Women's Medical college", 
        'ES': "St. Paul's Hospital Millennium Medical College", 
        'IN': 'National Institute of Cholera and Enteric Diseases', 
        'NK': 'Murtala Mohammed Specialist Hospital',
        'NN': 'National Hospital Abuja', 
        'NW': 'WUSE District Hospital',
        'PP': 'Pakistan Institute of Medical Sciences', 
        'PC': 'Bhara Kahu Health Centre',
        'RK': 'Kabgayi Hospital',
        'ZAT': 'Tygerberg Hospital',
        'RU': 'University Central Hospital of Kigali'
    }

    df['hospital_name'] = df['Site'].map(hosp_dictionary)
    hosp_code = pd.read_csv("FILEPATH")

    df = df.merge(hosp_code[['hospital_name', 'hospital_type']], how='left', on='hospital_name')

    assert df[['hospital_name', 'hospital_type']].notnull().values.all()

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

    mic_cols = [i for i in df.columns if 'MIC' in i]
    other_ast = ['Ampicillin2', 'Oxacillin2', 'Flucloxacillin2', 'Levofloxacin2', 'Ciprofloxacin2', 'Gentamicin2', 'Amikacin2', 'Tobramycin2', 
    'Tigecycline2','Minocycline2', 'Rifampicin2', 'Vancomycin2', 'Azithromycin2', 'Linezolide2']

    ast_cols = mic_cols+other_ast

    id_cols = [col for col in list(df) if col not in ast_cols]

    df = df.melt(id_vars=id_cols, value_vars=(mic_cols+other_ast))
    

    assert 'value' in df.columns

    df = df[~df.value.isna()]
    
    df['tmp'] = df.value.apply(extract_number)
    df['number'] = df['tmp'].apply(lambda x: x[0])
    df['number_val'] = df['tmp'].apply(lambda x: x[1])
    
    df['method'] = 'MIC'
    df['raw_antibiotic'] = df.variable.apply(lambda x: x[:-3] if x in mic_cols else x[:-1]) 
    
    return df
    
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
    

def create_antibiotic_resistance_interpretation(df):
    resistance_cols = []

    gpb_cols =['LVX2', 'CIP2', 'GEN2', 'AMK2',
                          'TOB2', 'TGC2', 'MIN2', 'RIF2', 'VAN2', 'AZM2',
                          'LZD2']
    gnb_cols = ['AMP', 'AMC', 'TZP', 'CRO', 'CTX', 'CAZ', 'FEP', 'IPM',
                'MEM', 'ETP', 'ATM', 'GEN', 'AMK', 'TOB', 'TGC', 'MIN',
                'FOF', 'LVX', 'CIP', 'CST']

    resistance_cols = gpb_cols + gnb_cols


    id_cols = [col for col in list(df) if col not in resistance_cols]

    df = df.melt(id_vars=id_cols, value_vars=resistance_cols,
                 var_name='raw_antibiotic', value_name='resistance')

    df.loc[df['resistance'].isin(['R', 'I']), 'resistance'] = 'resistant'
    df.loc[df['resistance'].isin(['0', 'S']), 'resistance'] = 'sensitive'


    df = df.loc[(df['resistance'].notnull()) & (df['resistance'] != 'ND'), ]

    return df



def clean_ab(df):
    
    df.raw_antibiotic = df.raw_antibiotic.str.lower()
    
    ab_dict = {'ampicilin':'ampicillin',
               'ampicilinmi':'ampicillin',
               'amoxicillin_clavulanicacid':'amoxicillin/clavulanic acid',
               'piperacilin_tazobactam':'piperacillin/tazobactam',
               'linezolide':'linezolid',
               'cefepimemi':'cefepime',
               'imipenemmi':'imipenem',
               'meropenemmi':'meropenem', 
               'ertapenemmi':'ertapenem', 
               'aztreonammi':'aztreonam', 
               'gentamicinmi':'gentamicin',          
               'amp':'ampicillin',
                'amc':'amoxicillin/clavulanic acid',
                'tzp':'piperacillin/tazobactam',
                'cro':'ceftriaxone', 
                'ctx':'cefotaxime',
                'caz':'ceftazidime',
                'fep':'cefepime',
                'ipm':'imipenem', 
                'mem':'meropenem',
                'etp':'ertapenem',
                'atm':'aztreonam',
                'gen':'gentamicin',
                'amk':'amikacin',
                'tob':'tobramycin',
                'tgc':'tigecycline',
                'min':'minocycline',
                 'fof':'fosfomycin',
                'lvx':'levofloxacin',
                'cip':'ciprofloxacin',
                'cst':'colistin',
               'lvx2':'levofloxacin',
            'cip2': 'ciprofloxacin',  
            'gen2':'gentamicin',
            'amk2': 'amikacin',
            'tob2':'tobramycin',
            'tgc2': 'tigecycline',
            'min2':'minocycline',
            'rif2':'rifampicin',
            'van2':'vancomycin',
            'azm2':'azithromycin',
            'lzd2':'linezolid'}

    df['raw_antibiotic'] = df.raw_antibiotic.apply(lambda x: ab_dict[x] if x in ab_dict.keys() else x)

    df = assign_breakpoints.validate_format_antibiotic_names(df)
    return df

def create_resistance_column_clsi_2023(df):

    df = clean_ab(df)  

    df['site_category'] = 'other'

    df = assign_breakpoints.assign_breakpoints(df, 'CLSI 2023')

    df['resistance'] = df.apply(lambda x: apply_resistance(x), axis=1)
    df.loc[df.resistance.isna(), 'resistance'] = 'unknown'

    return df

def clean_SOURCE():
    
    df = read_in_data()
    df = format_age_group_id(df)
    df = format_sex_id(df)
    df = format_location_id(df)
    df = mark_deaths(df)
    df = create_cause_cols(df)
    df = create_specimen_pathogen_col(df)

    df = create_hosp(df)
    df = add_hospital_type(df)
    df["year_id"] = pd.to_datetime(df["dateofadmissionintohospit"]).apply(lambda x: x.year)
    df.loc[df['year_id'].isna(), 'year_id'] = pd.to_datetime(df["startdatestudyparticipati"]).\
        apply(lambda x: x.year)
    df.loc[df['year_id'].isna(), 'year_id'] = pd.to_datetime(df["dateatinfantoutcome"]).\
        apply(lambda x: x.year)
    df.loc[df['year_id'].isna(), 'year_id'] = pd.to_datetime(df["babydob"]).\
        apply(lambda x: x.year)

    df = create_nid(df)

    df["cases"] = 1
    df['sample_id'] = 'SOURCE' + df['Anon_Isolate']
    df_I = create_antibiotic_resistance_interpretation(df)
    df = get_ast_values(df)
    
   
    df_I = clean_ab(df_I)

    df = create_resistance_column_clsi_2023(df)

    df_I.raw_pathogen = df_I.raw_pathogen.str.lower().str.replace(' ', '_')

    keep_cols =['sample_id', 'age_group_id', 'location_id', 'sex_id', 'year_id',
       'cases', 'raw_specimen', 'raw_pathogen',  'raw_antibiotic', 
                'nid', 'hosp', 'hospital_name', 'hospital_type', 'deaths']

    df_I = df_I[keep_cols+['resistance']]
    df_I.rename(columns={'resistance':'resistance_old'},inplace=True)

    df = df[keep_cols +['resistance','number','number_val', 'disk_dose','method', 'site_category']]

    m_cols = ['sample_id', 'age_group_id', 'location_id', 'sex_id', 'year_id',
           'cases', 'raw_specimen', 'raw_pathogen', 'raw_antibiotic','nid', 'hosp', 'hospital_name', 'hospital_type', 'deaths']

    df_I = df_I.drop_duplicates(m_cols)

    df = df.drop_duplicates(m_cols)

    df = df.merge(df_I, how='left', on=m_cols, validate='1:1')

    df['interpretation_method'] = df.apply(lambda x:'CLSI_2023' if x.resistance!='unknown' else 'other', axis=1)
    df['resistance']= df.apply(lambda x: x.resistance_old if x.resistance=='unknown' else x.resistance, axis=1)

    df.loc[df.resistance.isna(), 'resistance'] = 'unknown'
    print_log_message(f'{len(df[df.resistance!="unknown"])/len(df)*100} % of the dataset has resistance values assigned')
    
    return df

if __name__ == '__main__':
    df = clean_SOURCE()

    demo_cols = ['nid', 'location_id', 'year_id', 'age_group_id', 'sex_id']
    cause_cols = ['cause', 'multiple_cause_1', 'multiple_cause_2']
    biology = ['raw_specimen', 'raw_pathogen', 'raw_antibiotic', 'resistance','interpretation_method']
    other_values = ['hosp', 'hospital_name', 'hospital_type', 'sample_id', 'cases', 'deaths']
    lowercase_cols = ['raw_specimen', 'raw_pathogen', 'raw_antibiotic'] 
    
    resistance_testing = ['method', 'number', 'number_val', 'site_category'] 
    df = df[demo_cols + biology + other_values + resistance_testing]

    df['raw_pathogen']  = df.raw_pathogen.str.replace('_',' ')
    print(df.raw_pathogen.unique())

    for col in lowercase_cols:
        df[col] = df[col].str.lower().str.strip()
        df[col].notnull().values.all()

    df = df.drop_duplicates()
    df.to_csv("FILEPATH", index=False)
    df = df.drop(['method', 'number', 'number_val', 'site_category'], axis=1)
    
    save_amr_data(df, phase='unmapped', ldir=L_DIR)
