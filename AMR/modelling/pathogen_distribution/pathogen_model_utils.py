import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import sys
import getpass
import warnings
import itertools
import re
from cod_prep.utils import (
    report_duplicates, report_if_merge_fail
)
import itertools
from cod_prep.claude.claude_io import Configurator
from cod_prep.downloaders import (
    add_age_metadata, get_current_location_hierarchy,
    add_location_metadata, get_cod_ages, pretty_print,
    get_pop, add_population, getcache_age_aggregate_to_detail_map,
    get_country_level_location_id, get_ages,
    prep_age_aggregate_to_detail_map
)
from mcod_prep.utils.causes import get_infsyn_hierarchy, get_all_related_syndromes
from cod_prep.utils import print_log_message, warn_if_merge_fail


CONF = Configurator()


def add_country_location_id(df, **cache_kwargs):
    country_locs = get_country_level_location_id(
        df.location_id.unique().tolist(),
        get_current_location_hierarchy(
            location_set_version_id=CONF.get_id("location_set_version"),
            **cache_kwargs
        )
    )
    df = df.merge(country_locs, how='left', on='location_id', validate='many_to_one')
    report_if_merge_fail(df, 'country_location_id', 'location_id')
    return df


def add_model_ages(df, agg_age_group_ids, allow_fail=False):
    age_group_table = get_ages()
    good_age_group_ids = df.age_group_id.unique().tolist()
    age_detail_map = prep_age_aggregate_to_detail_map(age_group_table, good_age_group_ids)\
        .query("agg_age_group_id in @agg_age_group_ids")
    df = df.merge(
        age_detail_map[['age_group_id', 'agg_age_group_id']], how='left',
        on='age_group_id', validate='many_to_one'
    )
    if allow_fail:
        warn_if_merge_fail(df, 'agg_age_group_id', 'age_group_id')
    else:
        report_if_merge_fail(df, 'agg_age_group_id', 'age_group_id')
    return df


def aggregate_data(df, cols, agg_age_group_ids, id_cols,
                   value_cols=['cases']):
    can_aggregate = {'age_group_id', 'sex_id', 'hosp', 'location_id'}
    assert set(cols) <= can_aggregate,\
        f"Don't know how to aggregate {set(cols) - can_aggregate}"
    if 'age_group_id' in cols:
        df = add_model_ages(df, agg_age_group_ids, allow_fail=False)
        df['age_group_id'] = df['agg_age_group_id']
    if 'sex_id' in cols:
        df['sex_id'] = 3
    if 'hosp' in cols:
        df['hosp'] = 'all'
    if 'location_id' in cols:
        df = add_country_location_id(df, force_rerun=False, block_rerun=True)
        df['location_id'] = df['country_location_id']
    df = df.groupby(id_cols, as_index=False).agg(cases=('cases', 'sum'),
        redist_flag=('redist_flag', 'max')).reset_index()
    return df


class PathogenCFR(object):
    """Read and apply pathogen CFRs"""
    id_cols = ['location_id', 'year_id', 'age_group_id', 'sex_id', 'hosp']

    def __init__(self, infectious_syndrome, agg_age_group_ids, cfr_use=None,
                 cfr_ready=None):
        self.infectious_syndrome = infectious_syndrome
        self.agg_age_group_ids = agg_age_group_ids
        self.cfr_use = cfr_use or {}
        self.cfr_ready = cfr_ready

        self.conf = Configurator()
        self.cache_kwargs = {'force_rerun': True, 'block_rerun': False, 'cache_results': True}
        self.validate_inputs()

    def validate_inputs(self):
        """Validate inputs"""
        assert self.conf.config_type == 'amr'
        self.infsyn = get_infsyn_hierarchy(infsyn_set_version = 'gram2')
        assert self.infectious_syndrome in np.append(self.infsyn.infectious_syndrome.unique(), ['L1.5_respiratory_non_lower', 'L1.5_MMO_encephalitis'])
        ages = get_ages(**self.cache_kwargs)
        assert set(self.agg_age_group_ids) <= set(ages.age_group_id)

    def get_cfrs(self):
        """Get CFR models for pathogen/syndromes"""
        self.cfr_cols = self.id_cols + ['pathogen']
        value_col = 'predict'
        if self.cfr_ready:
            parent_dir = Path("FILEPATH")
            concatenated_df = pd.DataFrame()
            for cfrfile in parent_dir.iterdir():
                if cfrfile.suffix == '.csv' and cfrfile.name.startswith(re.sub(r'L\d+\.?\d*_', '', self.infectious_syndrome)):
                    df = pd.read_csv(f'FILEPATH/{cfrfile.name}')
                    concatenated_df = pd.concat([df, concatenated_df])
            df = concatenated_df
            assert df.notnull().values.all()
            if 'hosp' not in df:
                self.cfr_cols.remove('hosp')
            if (df.sex_id == 3).all():
                self.cfr_cols.remove('sex_id')

            if 'hosp' in df and 'unknown' not in df['hosp'].unique():
                warnings.warn(
                    "No hosp = 'unknown', averaging to get something"
                )
                df = df.append(
                    df.assign(hosp='unknown')
                    .groupby(self.cfr_cols, as_index=False)[value_col].mean()
                )
            if 'all' not in df.pathogen.unique():
                warnings.warn(
                    "No pathogen = 'all', assuming that other means all"
                )
                df = df.append(
                    df.query("pathogen == 'other'").assign(pathogen='all')
                )

            for pathogen, use_pathogen in self.cfr_use.items():
                df = df.loc[df.pathogen != pathogen]
                df = df.append(
                    df.query(f"pathogen == '{use_pathogen}'")
                    .assign(pathogen=pathogen)
                )

            df = df.drop_duplicates()
            report_duplicates(df, self.cfr_cols)
            df = df.rename(columns={value_col: 'cfr'})
            df = df[self.cfr_cols + ['cfr']]
            self.cfr = df.copy()
        else:
            self.cfr = pd.DataFrame(columns=self.cfr_cols + ['cfr'])

    def apply_cfrs(self, df, mark_rows=False):
        cfr = self.cfr.rename(columns={
            'location_id': 'country_location_id',
            'age_group_id': 'agg_age_group_id'
        })
        df = add_model_ages(df, self.agg_age_group_ids, allow_fail=False)
        df = add_country_location_id(df, **self.cache_kwargs)
        merge_cols = [
            c for c in self.cfr_cols if c not in ['location_id', 'age_group_id']
        ] + ['country_location_id', 'agg_age_group_id']
        df = df.merge(
            cfr[merge_cols + ['cfr']], how='left', validate='many_to_one',
            on=merge_cols
        )
        print_log_message(f"Filling {df.loc[df.cfr.isnull(), merge_cols].drop_duplicates()}")
        df = df.merge(
            cfr.loc[cfr.pathogen == 'all'].drop('pathogen', axis='columns')[
                [c for c in merge_cols if c != 'pathogen'] + ['cfr']
            ], how='left', validate='many_to_one',
            on=[c for c in merge_cols if c != 'pathogen']
        )
        df['cfr'] = df['cfr_x'].fillna(df['cfr_y'])
        if self.cfr_ready:
            report_if_merge_fail(df, 'cfr', merge_cols)
        else:
            df['cfr'] = df['cfr'].fillna(1)
        if mark_rows:
            df['using_cfrs'] = False
            df.loc[df.cases.isnull(), 'using_cfrs'] = True
        df['cases'] = df['cases'].fillna(df['deaths'] / df['cfr'])
        assert df['cases'].notnull().all()
        df = df.drop(
            ['agg_age_group_id', 'country_location_id', 'cfr',
             'cfr_x', 'cfr_y'], axis='columns'
        )
        return df
