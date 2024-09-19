from cod_prep.claude.configurator import Configurator
from pathlib import Path
from cod_prep.claude.claude_io import makedirs_safely
from cod_prep.utils import print_log_message, cod_timestamp
from amr_prep.utils.misc import get_prepped_csbg_universe
from amr_prep.utils.pathogen_redistributor import PathogenRedistributer
from amr_prep.utils.pathogen_sweeper import PathogenSweeper
import pandas as pd
import numpy as np
import glob
import os
from cod_prep.utils import wrap
import multiprocessing
from multiprocessing import Pool
from functools import partial
import tqdm
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

CONF = Configurator()


def read_results_wrapper(burden_type, cause_id, year_id, infectious_syndrome,
                         pathogen):
    return AmrResult(
        process='calculate_amr_burden',
        burden_type=burden_type,
        year_id=year_id,
        cause_id=cause_id,
        infectious_syndrome=infectious_syndrome,
        pathogen=pathogen
    ).read_results()


class AmrResult():
    burden_types = ["fatal", "nonfatal"]
    processes = {
        "split_sepsis_syndrome": {
            'valid_burden_types': ['fatal'],
            'valid_params': ['cause_id', 'infectious_syndrome']
        },
        "calculate_mi_ratios": {
            'valid_burden_types': ['nonfatal'],
            'valid_params': ['infectious_syndrome'],
        },
        "split_pathogen": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['cause_id', 'infectious_syndrome', 'pathogen']
        },
        "calculate_ylds_per_case": {
            'valid_burden_types': ['nonfatal'],
            'valid_params': ['infectious_syndrome']
        },
        "calculate_amr_props": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['pathogen']
        },
        "calculate_amr_burden": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['cause_id', 'infectious_syndrome', 'pathogen']
        },
        "aggregate_cause_syndrome": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['cause_id', 'infectious_syndrome', 'pathogen']
        },
        "summarize_burden": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['cause_id', 'infectious_syndrome', 'pathogen']
        },
    }

    conf = Configurator()

    def __init__(self, process: str, burden_type: str, year_id: int,
                 cause_id: int = None, infectious_syndrome: str = None,
                 pathogen: str = None, abx_class: str = None,
                 save_draws: bool = False, gram_version: str = None,
                 archive_version: str = None):
        self.process = process
        self.burden_type = burden_type
        self.year_id = year_id
        if (save_draws) & (process == 'summarize_burden'):
            self.draws_suffix = "_draws"
        else:
            self.draws_suffix = ""
        self.optional_params = {
            'cause_id': cause_id or "",
            'infectious_syndrome': infectious_syndrome or "",
            'pathogen': pathogen or "",
            'abx_class': abx_class or "",
        }
        assert self.conf.config_type == 'amr'
        self.gram_version = gram_version
        self.archive_version = archive_version
        self.construct_results_path()

    def construct_results_path(self):
        if self.process not in AmrResult.processes:
            raise NotImplementedError
        if self.burden_type not in AmrResult.burden_types:
            raise NotImplementedError

        assert self.burden_type in AmrResult.processes[self.process]['valid_burden_types']
        for param, val in self.optional_params.items():
            if param in AmrResult.processes[self.process]['valid_params']:
                assert val != "", f"You must pass a valid {param} for process {self.process}"
                if param == 'cause_id':
                    self.optional_params[param] = int(val)
            else:
                assert val == "", f"You cannot pass {param} for process {self.process}"

        if self.process == 'aggregate_cause_syndrome':
            subfolder = 'calculate_amr_burden'
        else:
            subfolder = self.process
        if self.archive_version is not None:
            subfolder = "FILEPATH"
        self.results_path = "FILEPATH"

    def return_filepath(self):
        return self.results_path

    def make_results_directory(self):
        makedirs_safely(str(self.results_path.parent))

    def clear_results(self):
        if self.archive_version is not None:
            raise AssertionError('Cannot clear an archive file')
        if self.results_path.exists():
            self.results_path.unlink()

    def write_results(self, df):
        if self.archive_version is not None:
            raise AssertionError('Cannot write to an archive file')
        self.make_results_directory()
        df.to_hdf(self.results_path, key='model', format='fixed', mode='w')

    def read_results(self):
        return pd.read_hdf(self.results_path, key='model')


def get_amr_results(gram_version, process, burden_type, year_id, cause_id=None,
                    infectious_syndrome=None, pathogen=None, abx_class=None, location_id=None,
                    age_group_id=None, sex_id=None, measure_id=None, metric_id=None, hosp=None,
                    counterfactual=None, exec_function=None, exec_function_args=None,
                    filter_continuously=True, draws=False, archive_version=None, keep_etpec=False):

    if gram_version == 'gram1':
        assert archive_version in ['RTR_M09D28', 'bug_only_w_flu_rsv'], \
            'If pulling gram1 results, must chose archive version'

    exec_function_args = exec_function_args or {}
    if process != 'summarize_burden':
        raise NotImplementedError
    assert burden_type in AmrResult.burden_types

    if (gram_version == 'gram2') & (burden_type == 'nonfatal') & (archive_version == 'v4'):
        csbd = get_prepped_csbg_universe(
            burden_type=burden_type,
            gram_version=gram_version,
            archive_version='v5',
            add_parents=True
        )
    else:
        csbd = get_prepped_csbg_universe(
            burden_type=burden_type,
            gram_version=gram_version,
            archive_version=archive_version,
            add_parents=True
        )
    if (burden_type == 'fatal') & (archive_version == 'v4') & ('L2_lower_respiratory_infection' == infectious_syndrome):
        infectious_syndrome = ['L2_lower_respiratory_infection', 'L2_mix_respiratory_infection']
    elif (infectious_syndrome == 'L2_mix_respiratory_infection'):
        AssertionError

    if not keep_etpec:
        csbd = csbd.loc[~csbd['pathogen'].isin(['enteropathogenic_escherichia_coli', 'enterotoxigenic_escherichia_coli']), ]

    csbd = csbd.loc[csbd.notnull().all(axis=1)]
    csbd_filters = {
        'cause_id': cause_id,
        'infectious_syndrome': infectious_syndrome,
        'pathogen': pathogen
    }
    for var, values in csbd_filters.items():
        if values is not None:
            csbd = csbd.loc[csbd[var].isin(wrap(values))]

    csbd = csbd[
        AmrResult.processes['summarize_burden']['valid_params']
    ].drop_duplicates()

    if (gram_version == 'gram1') and (archive_version == 'RTR_M09D28'):
        print('remove flu and rsv from gram1 pathogen only results')
        csbd = csbd.loc[
            ~(
                (csbd['infectious_syndrome'].isin(['respiratory_infectious','all'])) &
                (csbd['pathogen'].isin(['flu', 'rsv']))
            ), 
        ]

    print_log_message(f"Reading {len(csbd) * len(wrap(year_id))} files")

    filters = {
        'abx_class': abx_class,
        'age_group_id': age_group_id,
        'sex_id': sex_id,
        'location_id': location_id,
        'measure_id': measure_id,
        'metric_id': metric_id,
        'hosp': hosp,
        'counterfactual': counterfactual
    }

    dfs = []
    i = 0
    for year in wrap(year_id):
        for index, row in csbd.iterrows():
            if i % 50 == 0:
                print(f"Reading file {i}")
            df = AmrResult(
                process, burden_type,
                year, save_draws=draws,
                gram_version=gram_version,
                archive_version=archive_version, **row.to_dict()
            ).read_results()
            if exec_function is not None:
                df = exec_function(df, **exec_function_args)
            if filter_continuously:
                for var, values in filters.items():
                    if values is not None:
                        df = df.loc[df[var].isin(wrap(values))]
            dfs.append(df)
            i += 1
    df = pd.concat(dfs, sort=False)

    if not filter_continuously:
        for var, values in filters.items():
            if values is not None:
                df = df.loc[df[var].isin(wrap(values))]
    return df


def validate_amr_data(df, phase: str):
    rules = pd.read_csv("FILEPATH")\
        .query(f"phase == '{phase}'")

    allowed_columns = rules.column_name.unique()
    if phase == 'unmapped':
        mcause_cols = set([col for col in df if 'multiple_cause_' in col])
    else:
        mcause_cols = {}
    allowed_columns = list(set(allowed_columns).union(mcause_cols))
    assert df.columns.isin(allowed_columns).all(),\
        f"{set(df) - set(allowed_columns)} columns not allowed"

    required_columns = rules.query("required == 1").column_name.unique()
    assert set(required_columns) <= set(df),\
        f"{set(required_columns) - set(df)} columns are missing"

    no_nulls = rules.query("no_nulls == 1").column_name.unique()
    no_nulls = list(set(df).intersection(
        set(no_nulls).union(mcause_cols)
    ))
    for col in no_nulls:
        assert df[col].notnull().all(), f'{col} column has null values'

    data_types = rules.loc[rules.data_type.notnull()].set_index(
        'column_name')['data_type'].to_dict()
    check_cols = list(set(df).intersection(
        set(data_types).union(mcause_cols)
    ))
    for col in check_cols:
        if col in mcause_cols:
            assert df[col].dtype == str or df[col].dtype == object,\
            print("Failed Column: ", col)
        else:
            if data_types[col] in ('integer', 'numeric'):
                assert col in df.select_dtypes(include=['number']).columns,\
                print("Failed Column: ", col)
            elif data_types[col] == 'string':
                assert df[col].dtype == str or df[col].dtype == object,\
                print("Failed Column: ", col)

    allowed_vals = rules.loc[rules.allowed_values.notnull()].set_index(
        'column_name')['allowed_values'].to_dict()
    allowed_vals = {col: vals.split(',') for col, vals in allowed_vals.items()}
    check_cols = list(set(df).intersection(set(allowed_vals)))
    for col in check_cols:
        assert df[col].isin(allowed_vals[col]).all(),\
            f"{set(df[col]) - set(allowed_vals[col])} not allowed in {col}"


def save_amr_data(df: pd.DataFrame, phase: str, ldir: str):

    assert 'nid' in df, "You must have a real NID"
    assert df['nid'].notnull().values.all(), "You must have a real NID"
    assert np.issubdtype(df['nid'].dtype, np.number), "You must have a numeric NID"
    assert (df.nid != -1).all(), "You must have a real NID"

    validate_amr_data(df, phase=phase)

    if 'SOURCE' not in ldir:
        assert ("FILEPATH" in ldir) or (
            "FILEPATH" in ldir),\
            "You specified a directory not in the 'FILEPATH'"
    assert phase in ['unmapped', 'mapped']
    ldir = Path(ldir)
    assert ldir.exists()
    mapping_dir = "FILEPATH"
    archive_dir = "FILEPATH"
    makedirs_safely(str(mapping_dir))
    makedirs_safely(str(archive_dir))
    file_name = f"standardized_{phase}.csv"
    archive_file_name = f"standardized_{phase}_{cod_timestamp()}.csv"

    df.to_csv(mapping_dir / file_name, index=False)

    df.to_csv(archive_dir / archive_file_name, index=False)


def find_mapped_filepaths(target_sources_dict):

    target_source_path_dict = {}
    for key, value in target_sources_dict.items():
        for path in Path(value).rglob('*'):
            filepath = str(path)
            if "FILEPATH" in filepath:
                target_source_path_dict.update({key: path})

    return target_source_path_dict


def get_amr_data(sources=None, location_id=None, year_id=None, specimen=None,
                 pathogen=None, abx_class=None, infectious_syndrome=None, active=True,
                 redistribute_pathogens=True, aggregate_pathogens_only=False,
                 drop_unuseables=True, drop_remainder=True, add_cutoff=None):

    arguments = locals()

    nid_metadata = pd.read_csv("FILEPATH")

    # subset to active sources only if active = True but not active_flu, that is not used indefinitely
    if active:
        active_cols = [c for c in nid_metadata if ('active' in c) & ('active_flu' not in c)]
        nid_metadata = nid_metadata.loc[~(nid_metadata[active_cols] == 0).all(axis=1)]

    print_log_message('Getting the right filepaths')
    if sources is not None:
        if type(sources) is not list:
            sources = [sources]
        nid_metadata = nid_metadata.loc[nid_metadata['source'].isin(sources), ]
    wanted_source_dirs = nid_metadata.set_index('source')['file_path'].to_dict()

    args = [
        [source, filepath, redistribute_pathogens, aggregate_pathogens_only,
        drop_unuseables, drop_remainder, add_cutoff]
        for source, filepath in wanted_source_dirs.items()
    ]
    print_log_message('Reading in AMR data with parallelization')
    if redistribute_pathogens:
        print_log_message('Will redistribute unknown pathogen genuses')
        if drop_unuseables:
            print_log_message('Dropping unuseables pathogens such as "none", "other", "pathogen_unspecified"')
        if add_cutoff is not None:
            print_log_message('Limiting viable props to samples ' + str(add_cutoff) + ' or more within a prop_id')
        if drop_remainder:
            print_log_message('Also drop any remaining unspecified genera without valid prop merge')
    if aggregate_pathogens_only and not redistribute_pathogens:
        print_log_message('Aggregating pathogens only, no redistribution')
        if drop_unuseables:
            print_log_message('Dropping unuseables pathogens such as "none", "other", "pathogen_unspecified"')
    with Pool(processes=9) as pool:
        df_list = pool.map(parallelize_for_get_amr_data, args)
    df = pd.concat(df_list)

    del arguments['sources']
    del arguments['active']
    del arguments['redistribute_pathogens']
    del arguments['aggregate_pathogens_only']
    del arguments['drop_remainder']
    del arguments['drop_unuseables']
    del arguments['add_cutoff']
    for key, value in arguments.items():
        if value is not None:
            if type(value) is not list:
                value = [value]
            df = df.loc[df[key].isin(value), ]

    if 'sample_id' in df.columns:
        df['sample_id'] = df['sample_id'].str.replace(u'\xa0', u'')

    return df

def parallelize_for_get_amr_data(args):
    source, filepath, redistribute_pathogens, aggregate_pathogens_only, \
        drop_unuseables, drop_remainder, add_cutoff = args

    mappath = "FILEPATH"
    if os.path.isfile(mappath):
        df = pd.read_csv(mappath, low_memory=False)
        df['source'] = source

    drop_cols = []
    if redistribute_pathogens:
        redistributor = PathogenRedistributer(
            drop_remainder=drop_remainder,
            drop_unuseables=drop_unuseables,
            add_cutoff=add_cutoff
        )
        df = redistributor.get_computed_dataframe(df)
        drop_cols = [
            'iso3', 'age_group_years_start',
            'age_group_years_end', 'path_to_top_parent', '_merge',
            'aggregation_method', 'package_number'
        ]
    elif not redistribute_pathogens and aggregate_pathogens_only:
        sweeper = PathogenSweeper(sweep_type='map_aggregates', drop_unuseables=drop_unuseables)
        df = sweeper.get_computed_dataframe(df)
        drop_cols = ['aggregation_method', 'package_number']
        df.rename(columns={"pathogen": 'pre_aggregate_pathogen', 'final_list': 'pathogen'}, inplace=True)

    drop_cols = list(set(drop_cols).intersection(set(df.columns)))
    df.drop(columns=drop_cols, inplace=True)

    return df