
import numpy as np
import pandas as pd
import sys
from cod_prep.downloaders.locations import *
from cod_prep.downloaders import *
from mcod_prep.utils.nids import *
from amr_prep.utils.amr_io import *
from amr_prep.utils.pathogen_sweeper import PathogenSweeper
from mcod_prep.utils.mcause_io import get_mcause_data

WRITE = True
OUTDIR = "FILEPATH"
ISH = pd.read_csv("FILEPATH", encoding='latin1')
INITIAL_GROUPBY_COLS = ['location_id', 'year_id', 'age_group_id', 'infectious_syndrome', 'pathogen']


def get_melt_mcause_data():

    df = get_mcause_data(
            phase='format_map',
            is_active=True,
            sub_dirs='infectious_syndrome',
            verbose=True,
            force_rerun=True,
            block_rerun=False
        )

    df = df.join(df['pathogen_from_cause'].str.split(',', expand=True).add_prefix('pathogen_'))
    df.drop(columns=['pathogen_from_cause', 'cause_infectious_syndrome'], inplace=True)
    id_vars = [col for col in df.columns if 'pathogen_' not in col]
    df_pm = pd.melt(df, id_vars=id_vars, value_name='pathogen')
    df_pm = df_pm.loc[(df_pm['pathogen'].notnull()) & (df_pm['pathogen'] != 'none'), ]

    df = df_pm.rename(columns={'admissions': 'cases'})
    df.loc[df['cases'].isna(), 'cases'] = df['deaths']
    df['deaths'].fillna(0, inplace=True)

    df = df.groupby(INITIAL_GROUPBY_COLS, as_index=False)['deaths', 'cases'].sum()

    return df


def get_subset_amr_data():

    df = get_amr_data(active=True, redistribute_pathogens=False)

    nid = pd.read_csv("FILEPATH")
    df = df.merge(nid, how='left', on='source', validate='many_to_one')
    df = df.loc[df['active_3a'] == 1, ]

    if 'sample_id' in df.columns:
        nsd = df.loc[df['sample_id'].isna(), ]
        sd = df.loc[df['sample_id'].notnull(), ]
        sd = sd.drop_duplicates(subset=[col for col in sd.columns if col not in ['abx_class', 'resistance']])
    df = pd.concat([nsd, sd])

    df = df.groupby(INITIAL_GROUPBY_COLS, as_index=False)['cases', 'deaths'].sum()

    return df


def create_wb_locations(df):
    wblh = get_current_location_hierarchy(release_id=9, location_set_id=26, location_set_version_id=1115)
    wblh.rename(columns={'region_name': 'WB Income Level'}, inplace=True)
    lh = get_current_location_hierarchy()
    lh = lh.merge(wblh[['location_id', 'WB Income Level']], on='location_id', how='left')                    
    lh = lh.loc[lh['WB Income Level'].notnull(), ['iso3', 'WB Income Level']].drop_duplicates()

    df = add_location_metadata(df, 'iso3')
    df = df.merge(lh[['iso3', 'WB Income Level']], on='iso3', how='left', validate='many_to_one')

    print(df.groupby('WB Income Level', as_index=False)['cases'].sum())
    df.loc[df['WB Income Level'] != 'World Bank High Income', 'WB Income Level'] = 'World Bank Middle and Low Income'

    df_all = df.copy()
    df_all['WB Income Level'] = 'All Incomes'
    df = pd.concat([df, df_all])
    
    return df


def simple_global_year_bin(df, midpoint=2005):
    df = df.loc[df['year_id'] >= 1980, ]
    df_all = df.groupby(
        ['WB Income Level', 'age_group_id', 'infectious_syndrome', 'pathogen'],
        as_index=False
    )['cases', 'deaths'].sum()
    df_all['year_bin'] = '1980-2022'
    
    df.loc[df['year_id'].between(1980, midpoint-1), 'year_bin'] = '1980-' + str(midpoint-1)
    df.loc[df['year_id'].between(midpoint, 2022), 'year_bin'] = str(midpoint) + '-2022'    
    df = pd.concat([df, df_all])
    
    return df


def simple_age_binning(df):
    df = add_age_metadata(df, ['age_group_years_start', 'age_group_years_end'])
    df_neo = df.loc[df['age_group_years_end'] <= 0.07671233, ]
    df_u1 = df.loc[df['age_group_years_end'] <= 1, ]
    df_u5 = df.loc[df['age_group_years_end'] <= 5, ]
    df_p5 = df.loc[(df['age_group_years_start'] >= 0.07671233) & (df['age_group_years_end'] <= 5), ]
    df_5p = df.loc[(df['age_group_years_start'] >= 5), ]
    df_all = df.copy()

    df_neo['age_group_name'] = 'Neonatal'
    df_u1['age_group_name'] = '<1 year'
    df_u5['age_group_name'] = 'Under 5'
    df_p5['age_group_name'] = 'Post Neonatal to 5'
    df_5p['age_group_name'] = '5 plus'
    df_all['age_group_name'] = 'All Ages'
    
    df = pd.concat([df_neo, df_u1, df_u5, df_p5, df_5p, df_all])

    return df


def use_level_1_syndrome(df):
    lvl1_syn = ISH[['infectious_syndrome', 'path_to_top_parent']]
    lvl1_syn[['0', '1', '2', '3', '4']] = lvl1_syn['path_to_top_parent'].str.split(',', expand=True)
    df.loc[df['infectious_syndrome'] == 'lri_bordetella/b_parapertussis', 'infectious_syndrome'] = 'lri_bordetella_pertussis/parapertussis'
    df.loc[df['infectious_syndrome'] == 'meningitis_haemophilus', 'infectious_syndrome'] = 'meningitis_haemophilus_influenzae'
    df.loc[df['infectious_syndrome'] == 'bsi_salmonellae_ints_non_typhoidal', 'infectious_syndrome'] = 'bsi_salmonella_ints_non_typhoidal'

    df = df.merge(lvl1_syn, how='left', on='infectious_syndrome', indicator=True)

    df.drop(columns=['infectious_syndrome'], inplace=True)
    df.rename(columns={'1': 'infectious_syndrome'}, inplace=True)
    
    return df


def subset_pathogen_gen_props(df, drop_deaths=True):
    sweeper = PathogenSweeper(sweep_type='map_aggregates')
    df = sweeper.get_computed_dataframe(df)
    df.rename(columns={"pathogen": 'pre_aggregate_pathogen', 'final_list': 'pathogen'}, inplace=True)
    df = df.loc[(df['package_number'].notnull()) & (df['aggregation_method'] != 'redistributable'), ]
    df['genus'] = df['pathogen'].str.split('_', 1, expand=True)[0]

    df = df.groupby(
        ['WB Income Level', 'infectious_syndrome', 'genus', 'pathogen', 'age_group_name', 'year_bin'],
        as_index=False
    )['cases', 'deaths'].sum()
    
    df[['total_cases_genus', 'total_deaths_genus']] = \
        df.groupby(['WB Income Level', 'infectious_syndrome', 'genus', 'age_group_name', 'year_bin'],
                   as_index=False)['cases', 'deaths'].transform(sum)
    df['prop_cases'] = df['cases'] / df['total_cases_genus']
    df['prop_deaths'] = df['deaths'] / df['total_deaths_genus']

    df['# species in prop'] = df.groupby(
        ['WB Income Level', 'infectious_syndrome', 'genus', 'age_group_name', 'year_bin'],
        as_index=False
    )['pathogen'].transform('nunique')

    df['prop_id'] = df.groupby(
        ['WB Income Level', 'infectious_syndrome', 'genus', 'age_group_name', 'year_bin']
    ).ngroup().add(1)

    df['prop_sum'] = df.groupby('prop_id')['prop_cases'].transform('sum')
    assert np.isclose(df['prop_sum'], 1).all()
       
    df.rename(columns={'cases': 'species_cases', 'deaths': 'species_deaths'}, inplace=True)
    
    if drop_deaths:
        df.drop(columns=['species_deaths', 'total_deaths_genus', 'prop_deaths'], inplace=True)

    return df


def modify_prop_artificial(df):
    gi_corynb = df.loc[
        (df['infectious_syndrome'] == 'L1_gastrointestinal_infection') & 
        (df['genus'] == 'corynebacterium'),
    ]
    gi_corynb['infectious_syndrome'] = 'L1_peritoneal_and_intra_abdomen_infection'

    df = pd.concat([df, gi_corynb])

    return df


def generation_pathogen_redistribution_proportions():
    df = get_subset_amr_data()
    
    df = df.loc[df['cases'] != 0, ]
    
    df = create_wb_locations(df)
    df = simple_global_year_bin(df, midpoint=2005)
    df = simple_age_binning(df)
    df = use_level_1_syndrome(df)
    df = subset_pathogen_gen_props(df, drop_deaths=True)
    df = modify_prop_artificial(df)

    return df


if __name__ == '__main__':
    df = generation_pathogen_redistribution_proportions()
    if WRITE:
        df.to_csv("FILEPATH", index=False)

