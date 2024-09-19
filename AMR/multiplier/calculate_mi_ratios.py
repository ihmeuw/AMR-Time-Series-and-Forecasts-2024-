import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path
from get_draws.api import get_draws
from db_queries import get_outputs
from cod_prep.downloaders import (
    get_current_location_hierarchy,
    get_cod_ages,
    add_population,
    get_current_cause_hierarchy,
    getcache_age_aggregate_to_detail_map
)
from cod_prep.utils import (
    report_if_merge_fail, print_log_message,
    report_duplicates
)
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import AmrResult
from amr_prep.utils.misc import drop_unmodeled_asc
from multiply.split_pathogen import PathogenSplitter


CACHE_KWARGS = {
    'force_rerun': False,
    'block_rerun': True,
    'cache_results': False
}

CONF = Configurator('standard')
RELEASE_ID = CONF.get_id("release")
LSV_ID = CONF.get_id('location_set_version')
CSV_ID = CONF.get_id('computation_cause_set_version')
CS_ID = CONF.get_id('computation_cause_set')
POP_RUN = CONF.get_id('pop_run')

PROXY_CAUSES = {
    'L2_meningitis': ['meningitis'],
    'L2_encephalitis': ['encephalitis'],
    'L2_endocarditis': ['cvd_endo'],
    'L1.5_myocarditis_pericarditis_carditis': ['cvd_cmp_myocarditis'],
    'L2_lower_respiratory_infection': ['lri'],
    'L2_upper_respiratory_infection': ['uri'],
    'L2_typhoid_paratyphoid_ints': ['intest_typhoid', 'intest_paratyph', 'intest_ints'],
    'L2_diarrhea': ['diarrhea'],
    'L2_skin_infection': ['skin_bacterial', 'skin_cellulitis', 'skin_decubitus', 'skin_other'],
    'L2_urinary_tract_infection': ['urinary_nephritis'],
    'L2_sexually_transmitted_infection': ['std', 'std_chlamydia', 'std_gonnorhea', 'std_syphilis'],
    'L2_blood_stream_infection': ['maternal_sepsis', 'neonatal_sepsis'],
    'L2_tb': ['tb'],
    'L2_hepatitis': ['hepatitis', 'hepatitis_a', 'hepatitis_b', 'hepatitis_c', 'hepatitis_e'],
    'L2_other_parasitic_infection': ['_ntd'],
    'L2_unspecified_site_infection': ['infectious'],
    'L2_OAMOIS': ['infectious']
}

CFR_SYNS = [
    'L2_meningitis',
    'L2_encephalitis',
    'L2_myelitis_meningoencephalitis_other',
    'L2_lower_respiratory_infection',
    'L1_peritoneal_and_intra_abdomen_infection',
    'L2_skin_infection',
    'L2_oral_infection',
    'L2_eye_infection',
    'L2_genital_infection',
    'L1_bone_joint_infection',
    'L2_urinary_tract_infection',
    'L2_blood_stream_infection',
 ]

locs = get_current_location_hierarchy(
    location_set_version_id=LSV_ID, **CACHE_KWARGS
)
ALL_CNTRY = list(locs.loc[locs['level'] == 3, 'location_id'].unique())
AGES = get_cod_ages(**CACHE_KWARGS)['age_group_id'].unique().tolist()


def draws_wrapper(source, causes, year_id, measure_id, metric_id, version_id,
                  downsample, num_workers):
    draws = get_draws(
        gbd_id_type=["cause_id"],
        gbd_id=causes,
        source=source,
        year_id=year_id,
        measure_id=measure_id,
        metric_id=metric_id,
        version_id=version_id,
        release_id=RELEASE_ID,
        location_id=ALL_CNTRY,
        sex_id=[1, 2],
        age_group_id=AGES,
        downsample=downsample,
        n_draws=CONF.get_id('draw_num'),
        num_workers=num_workers
    )
    return draws


def get_mi_ratios(num_workers, year_id, infectious_syndrome):

    proxy_causes = PROXY_CAUSES[infectious_syndrome]
    ch = get_current_cause_hierarchy(
        cause_set_id=CS_ID, cause_set_version_id=CSV_ID,
        **CACHE_KWARGS
    ).set_index("acause")['cause_id'].to_dict()
    proxy_cause_ids = [ch[acause] for acause in proxy_causes]

    print_log_message(f"Pulling incidence for causes {proxy_cause_ids}")
    incd = draws_wrapper(
        source='como',
        causes=proxy_cause_ids,
        year_id=year_id,
        measure_id=6,
        metric_id=3,
        version_id=CONF.get_id("como_version"),
        downsample=True,
        num_workers=num_workers
    )

    print_log_message(f"Pulling deaths for causes {proxy_cause_ids}")
    mort = draws_wrapper(
        source='codcorrect',
        causes=proxy_cause_ids,
        year_id=year_id,
        measure_id=1,
        metric_id=1,
        version_id=CONF.get_id("codcorrect_version"),
        downsample=True,
        num_workers=num_workers
    )

    print_log_message("Calculating MI ratios...")
    demo_cols = ['year_id', 'location_id', 'age_group_id', 'sex_id']
    incd = add_population(incd, pop_run_id=POP_RUN, **CACHE_KWARGS)
    report_if_merge_fail(incd, 'population', demo_cols)

    draw_cols = [col for col in incd if 'draw_' in col]
    incd.loc[:, draw_cols] = incd.loc[:, draw_cols].multiply(
        incd.loc[:, 'population'], axis="index")

    incd = incd[demo_cols + ['cause_id'] + draw_cols]
    mort = mort[demo_cols + ['cause_id'] + draw_cols]

    pe = get_outputs(
        topic='cause',
        release_id=CONF.get_id('release'),
        compare_version_id=CONF.get_id('compare_version'),
        year_id=year_id,
        location_id=ALL_CNTRY, 
        age_group_id=AGES,
        measure_id=[1, 6],
        sex_id=[1, 2],
        metric_id=1,
        cause_id=proxy_cause_ids
    )
    pei = pe.loc[pe['measure_id'] == 6, demo_cols + ['cause_id', 'val']]
    ped = pe.loc[pe['measure_id'] == 1, demo_cols + ['cause_id', 'val']]
    pei.rename(columns={'val': 'pe_incd'}, inplace=True)
    ped.rename(columns={'val': 'pe_mort'}, inplace=True)
    incd = incd.merge(pei, on=demo_cols + ['cause_id'], how='left', validate='one_to_one')
    mort = mort.merge(ped, on=demo_cols + ['cause_id'], how='left', validate='one_to_one')

    age_meta_df = get_cod_ages(**CACHE_KWARGS)
    cause_meta_df = get_current_cause_hierarchy(
        cause_set_id=CS_ID, cause_set_version_id=CSV_ID,
        **CACHE_KWARGS
    )
    mort = drop_unmodeled_asc(mort, cause_meta_df, age_meta_df, 'fatal')
    incd = drop_unmodeled_asc(incd, cause_meta_df, age_meta_df, 'nonfatal')

    merge_cols = demo_cols + ['cause_id']
    mir = pd.merge(
        mort, incd, how='outer', on=merge_cols,
        suffixes=('_mortality', '_incidence'), indicator=True
    )
    mort_draws = [f'draw_{i}_mortality' for i in range(0, CONF.get_id('draw_num'))]
    incd_draws = [f'draw_{i}_incidence' for i in range(0, CONF.get_id('draw_num'))]

    assert (mir._merge != 'left_only').all()
    if 383 in mir['cause_id'].unique().tolist():
        mir = mir.loc[~((mir['cause_id'] == 383) & ~(mir['age_group_id'].isin([2, 3]))), ]
    mir['incidence_only'] = (mir._merge == 'right_only') | (mir[mort_draws] == 0).any(axis=1)
    assert mir.loc[mir.incidence_only, 'cause_id'].isin(
        [319, 320, 959, 297, 395, 396, 399, 400, 401, 402, 403, 404, 322]).all()
    assert mir.loc[mir._merge == 'both'].notnull().values.all()

    mir[mort_draws + ['pe_mort']] = mir[mort_draws + ['pe_mort']].fillna(0)
    assert mir.notnull().values.all()

    mir[draw_cols] = pd.DataFrame(
        mir.filter(regex='.*_incidence').to_numpy()
        / mir.filter(regex='.*_mortality').to_numpy(),
        index=mir.index
    )
    mir['mi_ratio'] = mir['pe_incd'] / mir['pe_mort']
    assert mir.loc[~mir.incidence_only, draw_cols].notnull().values.all()
    assert (mir.loc[~mir.incidence_only, draw_cols] != np.inf).values.all()
    mir = mir[
        demo_cols + ['cause_id', 'incidence_only'] + incd_draws + ['pe_mort', 'pe_incd', 'mi_ratio']
    ]
    report_duplicates(mir, demo_cols + ['cause_id'])
    return mir


def get_cfr(year_id, infectious_syndrome):
    print_log_message("Reading in CFRs...")
    base_dir = "FILEPATH"
    cfr = pd.read_csv("FILEPATH")
    cfr['syndrome'] = cfr['syndrome'].str.replace('_all', '')
    if infectious_syndrome in ['L2_myelitis_meningoencephalitis_other', 'L2_encephalitis']:
        cfr_mmo_enceph = cfr.loc[cfr['syndrome'] == 'MMO_encephalitis', ]
        cfr_mmo_enceph['infectious_syndrome'] = infectious_syndrome
        cfr = cfr.loc[cfr['syndrome'] != 'MMO_encephalitis', ]
    cfr_syn_rename = {
        'L2_meningitis': 'meningitis',
        'L2_encephalitis': 'MMO_encephalitis',
        'L2_myelitis_meningoencephalitis_other': 'MMO_encephalitis',
        'L2_lower_respiratory_infection_hosp': 'lower_respiratory_infection_hospital',
        'L2_skin_infection': 'skin_infection',
        'L2_oral_infection': 'oral_infection',
        'L2_eye_infection': 'eye_infection',
        'L2_urinary_tract_infection_hosp': 'urinary_tract_infection_hospital',
        'L2_genital_infection': 'genital_infection',
        'L2_blood_stream_infection': 'blood_stream_infection',
        'L1_bone_joint_infection': 'bone_joint_infection',
        'L1_peritoneal_and_intra_abdomen_infection': 'peritoneal_and_intra_abdomen_infection',
        'L2_lower_respiratory_infection_comm': 'lower_respiratory_infection_community',
        'L2_urinary_tract_infection_comm': 'urinary_tract_infection_community',
    }
    inv_cfr_syn_rename = {v: k for k, v in cfr_syn_rename.items()}
    cfr['infectious_syndrome'] = cfr['syndrome'].map(inv_cfr_syn_rename)
    if infectious_syndrome in ['L2_myelitis_meningoencephalitis_other', 'L2_encephalitis']:
        cfr = pd.concat([cfr, cfr_mmo_enceph])
    if infectious_syndrome == 'L2_lower_respiratory_infection':
        cfr = cfr.loc[cfr['infectious_syndrome'].isin(['L2_lower_respiratory_infection_comm', 'L2_lower_respiratory_infection_hosp']), ]
        cfr.loc[cfr['infectious_syndrome'] == 'L2_lower_respiratory_infection_comm', 'hosp'] = 'community'
        cfr.loc[cfr['infectious_syndrome'] == 'L2_lower_respiratory_infection_hosp', 'hosp'] = 'hospital'
    elif infectious_syndrome == 'L2_urinary_tract_infection':
        cfr = cfr.loc[cfr['infectious_syndrome'].isin(['L2_urinary_tract_infection_comm', 'L2_urinary_tract_infection_hosp']), ]
        cfr.loc[cfr['infectious_syndrome'] == 'L2_urinary_tract_infection_comm', 'hosp'] = 'community'
        cfr.loc[cfr['infectious_syndrome'] == 'L2_urinary_tract_infection_hosp', 'hosp'] = 'hospital'
    else:
        cfr = cfr.loc[cfr['infectious_syndrome'] == infectious_syndrome, ]
    if 'syndrome' in cfr.columns:
        cfr.drop(columns=['syndrome'], inplace=True)
    assert cfr.notnull().values.all()
    cfr = cfr.loc[cfr['year_id'] == year_id, ]

    if 'sex_id' not in cfr.columns:
        cfr['sex_id'] = 3
    cfr['pathogen'] = 'all'
    cfr.drop(columns=['haqi', 'age_name'], inplace=True)
    cfr.rename(columns={'predict': 'cfr'}, inplace=True)
    cfr['1/cfr'] = 1 / cfr['cfr']
    return cfr


def main(num_workers, year_id, infectious_syndrome):
    dfs = []
    if infectious_syndrome in CFR_SYNS:
        dfs.append(get_cfr(year_id, infectious_syndrome).assign(
            infectious_syndrome=infectious_syndrome
        ))
    elif infectious_syndrome in PROXY_CAUSES:
        dfs.append(
            get_mi_ratios(num_workers, year_id, infectious_syndrome).assign(
                infectious_syndrome=infectious_syndrome
            )
        )
    df = pd.concat(dfs, sort=False)
    AmrResult(
        process='calculate_mi_ratios',
        burden_type='nonfatal',
        year_id=year_id,
        infectious_syndrome=infectious_syndrome
    ).write_results(df)
    print_log_message("Done!")


if __name__ == '__main__':
    num_workers = int(sys.argv[1])
    task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        base_dir = "FILEPATH"
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        infectious_syndrome = str(task_row['infectious_syndrome'])
    else:
        year_id = int(sys.argv[2])
        infectious_syndrome = str(sys.argv[3])
    main(num_workers, year_id, infectious_syndrome)
