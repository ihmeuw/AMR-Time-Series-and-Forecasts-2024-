import os
import argparse
import glob
import pandas as pd
import warnings
from cod_prep.utils import print_log_message, clean_icd_codes
from cod_prep.downloaders import get_map_version
from cod_prep.claude.configurator import Configurator
from cod_prep.claude.claude_io import makedirs_safely
from mcod_prep.mcod_mapping import MCoDMapper
import numpy as np
from multiprocessing import Pool
from functools import partial
from amr_prep.utils.amr_io import save_amr_data
from mcod_prep.utils.misc import clean_col_strings

this_dir = "FILEPATH"
repo_dir = "FILEPATH"

MAP_DIR = "FILEPATH"
L_DIR = "FILEPATH"
CONF = Configurator()
DEFAULT_CAUSE = "none"
DEFAULT_RESISTANCE = "unknown"


def get_formatted_data(source_dir, arg):

    print_log_message('Getting unmapped data from:')
    print(source_dir)

    files = glob.glob("FILEPATH")
    if any('standardized_unmapped.csv' in file for file in files):
        if arg == 'new':
            if any('standardized_mapped.csv' not in file for file in files):
                df = pd.read_csv("FILEPATH")
        else:
            df = pd.read_csv("FILEPATH")
    else:
        return None
    df = df.rename(columns={'specimen': 'raw_specimen'})
    if 'sample_id' in df:
        assert not ((df.sample_id.dtype == int) or (df.sample_id.dtype == float)),\
            "Need to do some kind of source + id appending here"

    return df


def map_pathogens(df, filepath):
    print_log_message('Mapping pathogens...')

    pathogen_map = pd.read_csv("FILEPATH", encoding='latin1')
    pathogen_map['raw_pathogen'] = pathogen_map['raw_pathogen'].str.replace('\r', '')
    pathogen_map['raw_pathogen'] = pathogen_map['raw_pathogen'].str.lower().str.strip()
    pathogen_map = pathogen_map[['raw_pathogen',  'pathogen', 'pathogen_type']].drop_duplicates()


    df['raw_pathogen'] = df['raw_pathogen'].str.lower().str.strip()
    df = df.merge(pathogen_map, how='left', on='raw_pathogen', validate='many_to_one')
    grant_bugs = [
        'escherichia_coli', 'klebsiella_pneumoniae', 'mycobacterium_tuberculosis',
        'neisseria_gonorrheae', 'non_typhoidal_salmonellae', 'salmonella_typhi_paratyphi',
        'shigella_spp', 'staphylococcus_aureus', 'streptococcus_pneumoniae'
    ]
    if ((df['pathogen'].loc[df['pathogen'].notna()].isin(grant_bugs).all())
            & (df['pathogen'].isna().any())):
        warnings.warn(
            "Not all data mapped, but seems only have grant pathogens in "
            + str(filepath) + " Contact NAME"
        )
    if ((df['pathogen'].loc[df['pathogen'].notna()].isin(grant_bugs).all())
            & (df['pathogen'].notna().any())):
        warnings.warn("Only have grant pathogens in " + str(filepath) + " Contact NAME")
    return df


def map_antibiotics(df):
    print_log_message('Mapping antibiotics...')
    antibiotic_map = pd.read_csv("FILEPATH")
    antibiotic_map['raw_antibiotic'] = antibiotic_map['raw_antibiotic'].str.lower().str.strip()
    antibiotic_map = antibiotic_map[['raw_antibiotic',  'abx_class']].drop_duplicates()
    df['raw_antibiotic'] = df['raw_antibiotic'].str.lower().str.strip()
    df = df.merge(antibiotic_map, how='left', on='raw_antibiotic', validate='many_to_one')

    return df


def map_specimens(df):
    print_log_message('Mapping specimen...')
    specimen_map = pd.read_csv("FILEPATH", encoding='latin1')
    specimen_map = clean_col_strings(specimen_map, 'raw_specimen', deep_clean=False)
    specimen_map = specimen_map[['raw_specimen',  'specimen', 'sterile']].drop_duplicates()

    df = clean_col_strings(df, 'raw_specimen', deep_clean=False)
    df = df.merge(specimen_map, how='left', on='raw_specimen', validate='many_to_one')

    print_log_message('Mapping specimen to syndrome...')
    specimen_to_syndrome = pd.read_csv("FILEPATH").set_index(
        "specimen")['infectious_syndrome'].to_dict()
    df['specimen_syndrome'] = df['specimen'].map(specimen_to_syndrome)

    return df


def map_cause_info(df, nid_metadata, source):
    cause_cols = [c for c in df if 'cause' in c]
    has_cause = len(cause_cols) > 0
    has_ucause = 'cause' in cause_cols
    if has_cause:
        code_system_id = nid_metadata.query(
            f"source == '{source}'").code_system_id.unique()[0]
        infsyn_set_version = nid_metadata.query(
            f"source == '{source}'").infsyn_set_version.unique()[0]
        if not np.isnan(code_system_id):
            code_system_id = int(code_system_id)
            assert infsyn_set_version == 'v1', "This dataset is marked with non-null "\
                "code system but non-default infectious hierarchy - option not "\
                "implemented"
            code_map_version_id = get_map_version(
                code_system_id, gbd_round_id=CONF.get_id("gbd_round"))
            if code_system_id in [1, 6]:
                for col in cause_cols:
                    df[col] = clean_icd_codes(df[col], remove_decimal=True)
            mapper = MCoDMapper(
                'infectious_syndrome', drop_p2=True,
                code_system_id=code_system_id, code_map_version_id=code_map_version_id,
                infsyn_set_version='gram2'
            )
        else:
            path_to_custom_map = "FILEPATH"
            mapper = MCoDMapper(
                'infectious_syndrome', drop_p2=True,
                path_to_ucause_map=path_to_custom_map,
                path_to_int_cause_map=path_to_custom_map,
                infsyn_set_version='gram2'
            )
        df = mapper.get_computed_dataframe(df, map_underlying_cause=has_ucause)
        df = df.drop(cause_cols, axis='columns')

    if has_cause or 'specimen_syndrome' in df:
        df['final_syndrome'] = np.NaN
        if has_cause:
            df['final_syndrome'] = df['infectious_syndrome']
        if 'specimen_syndrome' in df:
            df.loc[
                (df.final_syndrome.isnull()
                 | df.final_syndrome.isin(MCoDMapper.remainder_syndromes)),
                'final_syndrome'
            ] = df['specimen_syndrome']
        if 'pathogen' in df:
            pathogen_to_syndrome = {
                'salmonella_typhi': 'L3_typhoid_fever',
                'salmonella_paratyphi': 'L3_paratyphoid_fever',
                'salmonella_typhi_paratyphi': 'L2_typhoid_paratyphoid_ints',
                'salmonella_ints_non_typhi/paratyphi': 'L3_ints',
                'mycobacterium_tuberculosis': 'L2_tb',
                'mycobacterium_leprae': 'L2_leprosy',
                'neisseria_meningitidis': 'L2_meningitis',
                'hepatitis_a_virus': 'L3_hepatitis_viral',
                'hepatitis_b_virus': 'L3_hepatitis_viral',
                'hepatitis_c_virus': 'L3_hepatitis_viral',
            }
            for pathogen, syndrome in pathogen_to_syndrome.items():
                df.loc[
                    (df.pathogen == pathogen)
                    & (df.final_syndrome.isnull()
                       | df.final_syndrome.isin(MCoDMapper.remainder_syndromes)),
                    'final_syndrome'
                ] = syndrome
        df['infectious_syndrome'] = df['final_syndrome']
        assert df['infectious_syndrome'].notnull().all()
    if 'cause_id' in df:
        assert df['cause_id'].notnull().all()

    df = syndrome_specimen_pathogen_overrides(df)

    df = df.drop(
        [
            c for c in df if c in
            ['specimen_syndrome', 'code_id', 'cause_infectious_syndrome',
             'final_syndrome']
        ],
        axis='columns'
    )

    if source == 'SOURCE':
        if 'all_possible_syndromes' in df.columns:
            df.drop(columns=['all_possible_syndromes'], inplace=True)

    if 'most_severe' in df.columns:
        df.drop(columns=['most_severe'], inplace=True)

    return df


def map_cause_info_in_parallel(df, num_workers, nid_metadata, source):
    print_log_message("Working on creating data chunks")
    samples = df['sample_id'].unique().tolist()
    groups = list(range(0, num_workers)) * (len(samples) // num_workers + 1)
    groups = groups[0:len(samples)]
    sample_to_group = dict(zip(samples, groups))
    df['group_id'] = df['sample_id'].map(sample_to_group)

    print_log_message("Checking for data with no cause info")
    cause_cols = [c for c in df if 'cause' in c]
    assert 'cause' in cause_cols
    assert len(cause_cols) > 0
    df['no_info'] = (df[cause_cols] == DEFAULT_CAUSE).all(axis=1)
    df['no_info'] = df.groupby('sample_id').no_info.transform(np.all)
    print(
        f"{df.no_info.sum()} out of {len(df)} rows have no cause info,"
        f" mapping separately..."
    )
    save_df = df.loc[df.no_info].copy().drop(cause_cols, axis='columns')
    save_df = map_cause_info(
        save_df, nid_metadata=nid_metadata,
        source=source
    )
    assert 'infectious_syndrome' in save_df
    save_df['cause_id'] = 919
    save_df['pathogen_from_cause'] = "none"
    df = df.loc[~df.no_info]

    data_chunks = [chunk for group, chunk in df.groupby('group_id')]
    _map_cause_info = partial(
        map_cause_info, nid_metadata=nid_metadata, source=source
    )

    print_log_message(f"Working on running cause mapping")
    pool = Pool(num_workers)
    df_list = pool.map(_map_cause_info, data_chunks)
    pool.close()
    pool.join()
    print_log_message("Done, putting chunks back together")
    df = pd.concat(df_list, axis='index', ignore_index=True, sort=True)
    df = df.append(save_df)
    df = df.drop(['group_id', 'no_info'], axis='columns')
    return df


def syndrome_specimen_pathogen_overrides(df):

    print_log_message('mixed respiratory overrides based on viral vs non-viral pathogen')
    not_viral = (df['pathogen_type'] != 'virus')
    mixed = (df['infectious_syndrome'] == 'L2_mix_respiratory_infection') | (df['infectious_syndrome'].str.contains('mri_'))
    df.loc[mixed & not_viral, 'infectious_syndrome'] = 'L2_lower_respiratory_infection'

    salmonella_remap = {
        'typhoid_salmonella_typhi': ['salmonella_typhi'],
        'paratyphoid_salmonella_para_typhi': ['salmonella_paratyphi'],
        'salmonellae_ints_non_typhoidal': ['salmonella_ints_non_typhi/paratyphi'],
        'L2_typhoid_paratyphoid_ints': ['salmonella_typhi', 'salmonella_ints_non_typhi/paratyphi', 'salmonella_paratyphi'],
        'L3_typhoid_fever': ['salmonella_typhi'],
        'L3_paratyphoid_fever': ['salmonella_paratyphi'],
        'L3_ints': ['salmonella_ints_non_typhi/paratyphi']
    }
    if ('specimen' in df.columns):
        print_log_message('Salmonella syndrome override if the pathogen is not salmonella related')
        for salmonella_syn, salmonella_pathogens in salmonella_remap.items():
            df.loc[
                (df['infectious_syndrome'] == salmonella_syn) &
                ~(df['pathogen'].isin(salmonella_pathogens)),
            'infectious_syndrome'] = df['specimen_syndrome']

    df.loc[df['pathogen'] == 'salmonella_typhi', 'infectious_syndrome'] = 'L3_typhoid_fever'
    df.loc[df['pathogen'] == 'salmonella_paratyphi', 'infectious_syndrome'] = 'L3_paratyphoid_fever'
    df.loc[df['pathogen'] == 'salmonella_ints_non_typhi/paratyphi', 'infectious_syndrome'] = 'L3_ints'

    return df


def map_all(arg, write=False):
    if write:
        print_log_message(
            "You passed write = True, this should only be done if your cleaning script"
            " has been approved in a PR review"
        )
    print_log_message("Searching for " + arg + " formatted files to map")
    nid_metadata = pd.read_csv("FILEPATH")

    if (arg != 'all') & (arg != 'new'):
        assert arg in nid_metadata['source'].unique()
        filepaths = nid_metadata.query(f"source == '{arg}'")['file_path'].unique().tolist()
    else:
        filepaths = nid_metadata['file_path'].unique().tolist()

    for filepath in filepaths:
        source = nid_metadata.query(f"file_path == '{filepath}'")['source'].unique()[0]
        df = get_formatted_data(filepath, arg)
        if df is None:
            print_log_message('No unmapped output found for: ' + filepath
                              + ' Has this been formatted?:')
            print_log_message(filepath)
            continue

        has_spec = False
        has_drug = False

        assert 'nid' in df, "You must have a real NID"
        assert np.issubdtype(df['nid'].dtype, np.number), "You must have a real NID"
        assert (df.nid != -1).all(), "You must have a real NID"

        if 'resistance' in df: 
            df.loc[~df["resistance"].isin(["resistant", "susceptible", "sensitive", "unknown"]), 'resistance'] = "unknown"
        if 'raw_pathogen' in df:
            df = map_pathogens(df, filepath)
        if 'raw_antibiotic' in df:
            df = map_antibiotics(df)
            has_drug = True
        if 'raw_specimen' in df:
            df = map_specimens(df)
            has_spec = True

        print_log_message(
            'Proceding to read out full mapped file, or unmapped specimen, drug, bugs')

        mapped_cols = {
            'specimen': 'raw_specimen',
            'abx_class': 'raw_antibiotic',
            'pathogen': 'raw_pathogen'
        }

        if not has_spec:
            del mapped_cols['specimen']
        if not has_drug:
            del mapped_cols['abx_class']

        any_unmapped = False

        for key, value in mapped_cols.items():
            unmapped = df[key].isna()
            if unmapped.values.any():
                print_log_message('Unmapped {} in the data, outputting a list of unmapped {}. '
                                  'Please update the map whenever possible'.format(value, value))
                unmapped_values = df.loc[unmapped, [value]].drop_duplicates()
                if len(unmapped_values)<5:
                    print_log_message(f'the following values were not mapped {unmapped_values}')
                any_unmapped = True
                if not os.path.exists("FILEPATH"):
                    makedirs_safely("FILEPATH")
                if os.path.exists("FILEPATH".format(value, source)):
                    os.remove("FILEPATH".format(value, source))
                unmapped_values.to_csv("FILEPATH".format(value, source), index=False)

        if source in ["SOURCES"]:
            df['cause'] = 'lower respiratory tract infections'

        if source == 'SOURCE':
            warnings.warn("This is SOURCE, make sure you're running with 50 threads")
            df = map_cause_info_in_parallel(df, 50, nid_metadata, source)
        else:
            df = map_cause_info(df, nid_metadata, source)

        if "hospital_type" in df.columns.tolist() or "hospital_name" in df.columns.tolist():
            assert "hospital_name" in df.columns.tolist()
            assert df[["hospital_type", "hospital_name"]].notnull(
            ).values.all(), "Not every hospital has been classified"

            df.drop("hospital_name", axis=1, inplace=True)
            df.loc[df['hospital_type'].isin(["unknown", "mixed"]), "hospital_type"] = "m/u"

            assert df.hospital_type.isin(["tertiary", "non-tertiary", "m/u"]).values.all(
            ), "there are more hospital classifications in the data than allowed"

        if not any_unmapped:
            print_log_message('All specimen, raw_pathogen, and raw_antibiotic mapped!')
            drop_cols = list(mapped_cols.values())
            if 'raw_antibiotic' in drop_cols:
                drop_cols.remove('raw_antibiotic')
            df = df.drop(columns=drop_cols)
            if write:
                print_log_message('Saving standardized_mapped file')
                save_amr_data(df, phase='mapped', ldir=filepath)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="remap all or only new data")
    parser.add_argument(
        "type", help="'all' to remap everything, 'new' for any datasets not mapped yet, "
        "or an individual source name"
    )
    write_parser = parser.add_mutually_exclusive_group(required=True)
    write_parser.add_argument(
        '--write', dest='write', action='store_true',
        help='Write the output file from mapping '
        '(only use once your cleaning script has been approved)'
    )
    write_parser.add_argument(
        '--no-write', dest='write', action='store_false',
        help='Do not write the output file from mapping (used for testing purposes)'
    )
    arg = parser.parse_args()

    map_all(arg.type, write=arg.write)
