import sys
import os
import numpy as np
import pandas as pd
from get_draws.api import get_draws
from db_queries import get_outputs
from cod_prep.downloaders import (
    get_current_location_hierarchy, get_cod_ages,
    get_current_cause_hierarchy, add_age_metadata
)
from cod_prep.utils import (
    print_log_message,
    report_duplicates
)
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import AmrResult
from amr_prep.utils.misc import drop_unmodeled_asc


CACHE_KWARGS = {
    'force_rerun': False,
    'block_rerun': True,
    'cache_results': False
}

CONF = Configurator('standard')
RELEASE_ID = CONF.get_id("release")
DECOMP_STEP = CONF.get_id("decomp_step")
LSV_ID = CONF.get_id('location_set_version')
CSV_ID = CONF.get_id('computation_cause_set_version')
CS_ID = CONF.get_id('computation_cause_set')
POP_RUN = CONF.get_id('pop_run')

PROXY_CAUSES = {
    'L2_meningitis': ['meningitis'],
    'L2_encephalitis': ['encephalitis'],
    'L2_myelitis_meningoencephalitis_other': ['encephalitis'],
    'L2_endocarditis': ['cvd_endo'],
    'L1.5_myocarditis_pericarditis_carditis': ['cvd_cmp_myocarditis'],
    'L2_lower_respiratory_infection': ['lri'],
    'L2_mix_respiratory_infection': ['lri'],
    'L2_upper_respiratory_infection': ['uri'],
    'L2_skin_infection': ['skin_infect', 'skin_cellulitis', 'skin_bacterial', 'skin_decubitus'],
    'L2_oral_infection': ['otitis'],
    'L2_eye_infection': ['otitis'],
    'L2_typhoid_paratyphoid_ints': ['intest_typhoid', 'intest_paratyph', 'intest_ints'],
    'L2_diarrhea': ['diarrhea'],
    'L2_sexually_transmitted_infection': ['std', 'std_chlamydia', 'std_gonnorhea', 'std_syphilis'],
    'L2_urinary_tract_infection': ['urinary_nephritis'],
    'L2_genital_infection': ['urinary_nephritis'],
    'L2_blood_stream_infection': ['maternal_sepsis', 'neonatal_sepsis'],
    'L2_tb': ['tb'],
    'L1_peritoneal_and_intra_abdomen_infection': ['digest_ileus'],
    'L1_bone_joint_infection': ['skin_infect', 'msk'],
    'L2_hepatitis': ['hepatitis', 'hepatitis_a', 'hepatitis_b', 'hepatitis_c', 'hepatitis_e'],
    'L2_other_parasitic_infection': ['_ntd'],
    'L2_unspecified_site_infection': ['infectious'],
    'L2_OAMOIS': ['infectious']
}
LH = get_current_location_hierarchy(
    location_set_version_id=LSV_ID, **CACHE_KWARGS
)
COUNTRIES = list(LH.loc[LH['level'] == 3, 'location_id'].unique())
AGES = get_cod_ages( **CACHE_KWARGS)['age_group_id'].unique().tolist()

def draws_wrapper(causes, year_id, ages, measure_id, downsample, num_workers):
    draws = get_draws(
        gbd_id_type=["cause_id"],
        gbd_id=causes,
        source='como',
        year_id=year_id,
        measure_id=measure_id,
        metric_id=3,
        release_id=RELEASE_ID,
        version_id=CONF.get_id("como_version"),
        location_id=COUNTRIES,
        sex_id=[1, 2],
        age_group_id=ages,
        num_workers=num_workers,
        downsample=downsample,
        n_draws=CONF.get_id("draw_num"),
    )
    return draws


def get_ylds_per_case(num_workers, year_id, infectious_syndrome):

    proxy_causes = PROXY_CAUSES[infectious_syndrome]
    ch = get_current_cause_hierarchy(
        cause_set_id=CS_ID, cause_set_version_id=CSV_ID,
        **CACHE_KWARGS
    ).set_index("acause")['cause_id'].to_dict()
    proxy_cause_ids = [ch[acause] for acause in proxy_causes]

    if 368 in proxy_cause_ids:
        ages = AGES + [22]
    else:
        ages = AGES

    print_log_message(f"Pulling incidence for causes {proxy_cause_ids}")
    incd = draws_wrapper(
        causes=proxy_cause_ids,
        year_id=year_id,
        ages=ages,
        measure_id=6,
        downsample=True,
        num_workers=num_workers
    )

    print_log_message(f"Pulling YLDs for causes {proxy_cause_ids}")
    ylds = draws_wrapper(
        causes=proxy_cause_ids,
        year_id=year_id,
        ages=ages,
        measure_id=3,
        downsample=True,
        num_workers=num_workers
    )
    demo_cols = ['year_id', 'location_id', 'age_group_id', 'sex_id']
    pe = get_outputs(
        topic='cause',
        release_id=CONF.get_id('release'),
        compare_version_id=CONF.get_id('compare_version'),
        year_id=year_id,
        location_id=COUNTRIES, 
        age_group_id=ages,
        sex_id=[1, 2],
        measure_id=[3, 6], 
        metric_id=1,
        cause_id=proxy_cause_ids
    )
    pei = pe.loc[pe['measure_id'] == 6, demo_cols + ['cause_id', 'val']]
    pey = pe.loc[pe['measure_id'] == 3, demo_cols + ['cause_id', 'val']]
    pei.rename(columns={'val': 'pe_incd'}, inplace=True)
    pey.rename(columns={'val': 'pe_ylds'}, inplace=True)
    incd = incd.merge(pei, on=demo_cols + ['cause_id'], how='left', validate='one_to_one')
    ylds = ylds.merge(pey, on=demo_cols + ['cause_id'], how='left', validate='one_to_one')

    print_log_message("Calculating YLDs per case...")
    demo_cols = ['year_id', 'location_id', 'age_group_id', 'sex_id']
    draw_cols = [col for col in incd if 'draw_' in col]
    incd = incd[demo_cols + ['cause_id'] + draw_cols + ['pe_incd']]
    ylds = ylds[demo_cols + ['cause_id'] + draw_cols + ['pe_ylds']]

    age_meta_df = get_cod_ages(**CACHE_KWARGS)
    cause_meta_df = get_current_cause_hierarchy(
        cause_set_id=CS_ID, cause_set_version_id=CSV_ID,
        **CACHE_KWARGS
    )
    incd = drop_unmodeled_asc(incd, cause_meta_df, age_meta_df, 'nonfatal')
    ylds = drop_unmodeled_asc(ylds, cause_meta_df, age_meta_df, 'nonfatal')

    merge_cols = demo_cols + ['cause_id']
    df = pd.merge(
        ylds, incd, how='outer', on=merge_cols, validate='one_to_one',
        suffixes=('_yld', '_incd'), indicator=True
    )
    assert (df._merge == 'both').all()

    if infectious_syndrome == 'L1_bone_joint_infection':
        df = add_age_metadata(df, 'age_group_years_end')
        bone = df.loc[(df['age_group_years_end'] > 5) & (df['cause_id'] == 626), ]
        skin_proxy = df.loc[(df['age_group_years_end'] <= 5) & (df['cause_id'] == 980), ]
        df = pd.concat([bone, skin_proxy])

    yld_draws = df.filter(regex='.*_yld').to_numpy()
    incd_draws = df.filter(regex='.*_incd').to_numpy()
    df[draw_cols + ['pe_ylds_case']] = pd.DataFrame(
        np.true_divide(
            yld_draws, incd_draws, out=np.zeros_like(yld_draws),
            where=incd_draws != 0
        ), index=df.index
    )
    assert df.notnull().values.all()
    assert (df != np.inf).values.all()
    df = df[demo_cols + ['cause_id'] + draw_cols + ['pe_ylds_case']]
    if infectious_syndrome == 'L2_blood_stream_infection':
        age_limits = tuple(
            cause_meta_df.query(f"cause_id == 368").iloc[0][
                ['yld_age_start', 'yld_age_end']
            ])
        add_ages = age_meta_df.query(
            f"simple_age < {age_limits[0]} | simple_age > {age_limits[1]}"
        ).age_group_id.unique().tolist()

        copy_df = df.query("age_group_id == 22 & cause_id == 368")
        df = df.append(pd.concat(
            [copy_df.assign(age_group_id=age) for age in add_ages]
        ))
        df = df.append(
            df.query("cause_id == 368 & sex_id == 2")
            .assign(sex_id=1)
        )

    report_duplicates(df, demo_cols + ['cause_id'])
    assert df.notnull().values.all()

    AmrResult(
        process='calculate_ylds_per_case',
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
    get_ylds_per_case(num_workers, year_id, infectious_syndrome)
