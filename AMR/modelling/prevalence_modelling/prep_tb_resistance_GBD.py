'''
leverage the GBD estimates of TB parent, TB-mdr, TB-xdr to calculate a resistance fraction 
'''

import numpy as np
import pandas as pd
from get_draws.api import get_draws
from db_queries import get_outputs
import sys
from cod_prep.downloaders import *
from cod_prep.utils import *
from cod_prep.claude.configurator import Configurator

CONF = Configurator('standard')

lh = get_current_location_hierarchy()
locs = lh.loc[lh['level'] == 3, 'location_id'].unique().tolist()
ages = [22]
sex = [3]
cause_ids = [297, 946, 947]
draw_cols = ["draw_" + str(x) for x in range(0, 100)]
dem_cols = ['year_id', 'location_id', 'age_group_id', 'sex_id']
out_dir = 'FILEPATH'


def prep_prevalence_resistance_tb_GBD(year):
    print_log_message('Pulling draws and point estimates')
    df_deaths, df_incidence = pull_GBD_numbers(year)

    print_log_message('Calculating fractions and writing out..')
    df_deaths = calculate_fraction(df_deaths)
    df_incidence = calculate_fraction(df_incidence)

    df_deaths['measure_id'] = 1
    df_incidence['measure_id'] = 6

    df = pd.concat([df_deaths, df_incidence])
    df.to_csv(out_dir + str(year) + '.csv', index=False)


def calculate_fraction(df):
    tbp = df.loc[df['cause_id'] == 297, dem_cols + draw_cols + ['val']].set_index(dem_cols)
    mdr = df.loc[df['cause_id'] == 946, dem_cols + draw_cols + ['val']].set_index(dem_cols)
    xdr = df.loc[df['cause_id'] == 947, dem_cols + draw_cols + ['val']].set_index(dem_cols)

    mdr = (mdr/tbp).reset_index()
    xdr = (xdr/tbp).reset_index()

    mdr['abx_class'] = 'mdr'
    xdr['abx_class'] = 'xdr'

    df = pd.concat([mdr, xdr])
    df.rename(columns={'val': 'point_estimate'}, inplace=True)

    assert df.notnull().values.all()
    assert (df[draw_cols + ['point_estimate']] >= 0).values.all()
    assert (df[draw_cols + ['point_estimate']] <= 1).values.all()

    return df


def pull_GBD_numbers(year):
    draws_deaths = get_draws(
        "cause_id",
        cause_ids,
        year_id=year,
        age_group_id=ages,
        sex_id=sex,
        source='codcorrect',
        measure_id=[1],
        metric_id=[1],
        location_id=locs,
        release_id=CONF.get_id('release'),
        version_id=CONF.get_id('codcorrect_version'),
        downsample=True,
        n_draws=CONF.get_id('draw_num'),
        num_workers=5
    )

    draws_incidence = get_draws(
        "cause_id",
        cause_ids,
        year_id=year,
        age_group_id=ages,
        sex_id=sex,
        source='como',
        measure_id=[6],
        metric_id=[3],
        location_id=locs,
        release_id=CONF.get_id('release'),
        version_id=CONF.get_id('como_version'),
        downsample=True,
        n_draws=CONF.get_id('draw_num'),
        num_workers=5
    )
    pop = get_population(
        release_id=CONF.get_id('release'),
        year_id=year,
        location_id=locs,
        age_group_id=ages,
        sex_id=sex
    )
    draws_incidence = draws_incidence.merge(
        pop,
        on=dem_cols,
        how='left',
        validate='many_to_one'
    )
    draws_incidence[draw_cols] = draws_incidence[draw_cols].mul(draws_incidence['population'], axis=0)

    pe = get_outputs(
        topic='cause',
        release_id=CONF.get_id('release'),
        compare_version_id=CONF.get_id('compare_version'),
        year_id=year,
        location_id=locs,
        age_group_id=ages,
        sex_id=sex,
        measure_id=[1, 6],
        metric_id=1, 
        cause_id=cause_ids
    )
    pe_deaths = pe.loc[pe['measure_id'] == 1, dem_cols + ['cause_id', 'val']]
    pe_incidence = pe.loc[pe['measure_id'] == 6, dem_cols + ['cause_id', 'val']]
    del pe

    df_deaths = draws_deaths.merge(pe_deaths, on=dem_cols + ['cause_id'], how='outer', validate='one_to_one') 
    df_incidence = draws_incidence.merge(pe_incidence, on=dem_cols + ['cause_id'], how='outer', validate='one_to_one') 
    df_deaths = df_deaths[dem_cols + ['cause_id'] + draw_cols + ['val']]
    df_incidence = df_incidence[dem_cols + ['cause_id'] + draw_cols + ['val']]
    
    assert ~df_deaths[dem_cols + ['cause_id']].duplicated().any()
    assert ~df_incidence[dem_cols + ['cause_id']].duplicated().any()

    del draws_deaths
    del draws_incidence
    del pe_deaths
    del pe_incidence

    return df_deaths, df_incidence


if __name__ == "__main__":
    year = int(sys.argv[1])
    prep_prevalence_resistance_tb_GBD(year)
