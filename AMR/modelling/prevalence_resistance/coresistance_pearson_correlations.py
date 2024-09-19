import pandas as pd
import numpy as np
import sys
import os
import argparse
import getpass
user = getpass.getuser()
amrpath = 'FILEPATH'.format(user)
if not amrpath in sys.path:
    sys.path.append('FILEPATH'.format(user))
from amr_prep.utils.amr_io import get_amr_data
from datetime import datetime

OUT_DIR = "FILEPATH"

amr_dat = get_amr_data(redistribute_pathogens = False)
for s in ['citrobacter','serratia','providencia','proteus','enterobacter','moraxella','morganella','shigella']:
    amr_dat.loc[amr_dat.pathogen.str.contains(s), 'pathogen'] = s + '_spp'

for s in ['campylobacter','listeria']:
    amr_dat.loc[amr_dat.pathogen.str.contains(s), 'pathogen'] = s

for s in ['group_b','group_a']:
    amr_dat.loc[amr_dat.pathogen.str.contains(s), 'pathogen'] = s + '_strep'
    
amr_dat.loc[amr_dat.pathogen.str.contains('enterococcus') & 
            ~amr_dat.pathogen.isin(['enterococcus_faecium','enterococcus_faecalis']), 'pathogen'] = 'enterococcus_spp'

amr_dat.loc[amr_dat.pathogen.str.contains('gonorrheae'), 'pathogen'] = 'neisseria_gonorrhoeae'
amr_dat.loc[amr_dat.pathogen.str.contains('baumanii'), 'pathogen'] = 'acinetobacter_baumannii'

resdat = amr_dat.loc[(amr_dat['resistance'].notna()) & (
    amr_dat['resistance'].isin(['sensitive', 'resistant', 'susceptible'])), :]
resdat['resistance'] = resdat['resistance'].str.replace('sensitive', 'susceptible')

del amr_dat

dupsamps = resdat.loc[resdat['sample_id'].duplicated(), 'sample_id'].unique()

multires = resdat.loc[resdat['sample_id'].isin(dupsamps), :]

studycombos = pd.read_csv(f'FILEPATH',encoding='latin-1')

multires = multires.merge(studycombos.loc[:, ['pathogen', 'abx_class']], on=[
                          'pathogen', 'abx_class'])
multires = multires.loc[:, ['pathogen', 'abx_class',
                            'cases', 'resistance', 'sample_id']].drop_duplicates()


def calc_class_resistance(df, groupcols):
    '''Fold multiple resistance values within a patient for
    the same antibiotic class into single value'''
    if df.groupby(groupcols).size().max() > 1:
        allcols = [col for col in df if col != 'resistance']
        computed = df.fillna('unknown').groupby(allcols, group_keys=False, as_index=False).apply(
            lambda x: pd.Series(
                {'resistance': 'resistant' if any(x['resistance'] == 'resistant')
                 else 'susceptible'}
            )
        ).reset_index()
        computed = computed.drop_duplicates()
        return(computed)
    else:
        return(df)


multires.loc[multires['resistance'] == 'susceptible', 'resistance'] = 0
multires.loc[multires['resistance'] == 'resistant', 'resistance'] = 1

multires = multires.loc[~multires['pathogen'].isin(
    ['neisseria_gonorrheae', 'shigella_spp', 'salmonella_paratyphi', 'group_a_strep',
     'salmonella_typhi', 'non_typhoidal_salmonellae']), :]

pathogen_list = []
abx_a_list = []
abx_b_list = []
nres_a_only = []
nres_b_only = []
nres_both = []
nsuscept_both = []

for bug in multires['pathogen'].unique():
    drugs = list(studycombos.loc[studycombos['pathogen'] == bug, 'abx_class'])
    for abxa in drugs:
        for abxb in [drug2 for drug2 in drugs if drug2 not in abxa]:
            print(datetime.now().strftime("%H:%M:%S") + bug + '-' + abxa + '-' + abxb)
            subset = multires.loc[(multires['pathogen'] == bug) & (
                multires['abx_class'].isin([abxa, abxb])), :]
            pivot = subset.pivot_table(
                index=['sample_id', 'pathogen', 'cases'], columns='abx_class',
                values='resistance').reset_index()
            pivot.rename(columns={abxa: 'abx_a', abxb: 'abx_b'}, inplace=True)

            pathogen_list.append(bug)
            abx_a_list.append(abxa)
            abx_b_list.append(abxb)
            nres_a_only.append(
                sum(pivot.loc[(pivot['abx_a'] > 0) & (pivot['abx_b'] == 0), 'cases']))
            nres_b_only.append(
                sum(pivot.loc[(pivot['abx_b'] > 0) & (pivot['abx_a'] == 0), 'cases']))
            nres_both.append(sum(pivot.loc[(pivot['abx_a'] > 0) & (pivot['abx_b'] > 0), 'cases']))
            nsuscept_both.append(
                sum(pivot.loc[(pivot['abx_a'] == 0) & (pivot['abx_b'] == 0), 'cases']))

coresist = pd.DataFrame({'pathogen': pathogen_list, 'abx_a': abx_a_list, 'abx_b': abx_b_list,
                         'n_res_a_only': nres_a_only, 'n_res_b_only': nres_b_only,
                         'n_res_ab': nres_both,
                         'n_suscept_ab': nsuscept_both})

coresist['tetra_corr'] = np.cos(np.pi / (1 + np.sqrt(coresist['n_suscept_ab']
                                                     * coresist['n_res_ab']
                                                     / coresist['n_res_a_only']
                                                     / coresist['n_res_b_only'])))

coresist.to_csv(f'{OUT_DIR}/tetrachoric_correlations.csv', index=False)

corrs = pd.read_csv(f'{OUT_DIR}/tetrachoric_correlations.csv')

corrs['prop_ab'] = corrs['n_res_ab'] / \
    (corrs['n_res_a_only'] + corrs['n_res_b_only'] + corrs['n_res_ab'] + corrs['n_suscept_ab'])
corrs['prop_a'] = (corrs['n_res_a_only'] + corrs['n_res_ab']) / (
    corrs['n_res_a_only'] + corrs['n_res_b_only'] + corrs['n_res_ab'] + corrs['n_suscept_ab'])
corrs['prop_b'] = (corrs['n_res_b_only'] + corrs['n_res_ab']) / (
    corrs['n_res_a_only'] + corrs['n_res_b_only'] + corrs['n_res_ab'] + corrs['n_suscept_ab'])

corrs['pearson'] = (corrs['prop_ab'] - corrs['prop_a'] * corrs['prop_b']) / (
    np.sqrt((corrs['prop_a'] * (1 - corrs['prop_a']) * corrs['prop_b'] * (1 - corrs['prop_b']))))

dropcombos = studycombos.loc[studycombos['kevin_catrin_exclude'] == 1, ['pathogen', 'abx_class']]

corrslimit = corrs.merge(dropcombos, left_on=['pathogen', 'abx_a'], right_on=[
                         'pathogen', 'abx_class'], how='outer', indicator=True)
corrslimit = corrslimit.loc[
    corrslimit['_merge'] == 'left_only', :].drop(columns=['abx_class', '_merge'])
corrslimit = corrslimit.merge(dropcombos, left_on=['pathogen', 'abx_b'], right_on=[
                              'pathogen', 'abx_class'], how='outer', indicator=True)
corrslimit = corrslimit.loc[
    corrslimit['_merge'] == 'left_only', :].drop(columns=['abx_class', '_merge'])
corrslimit = corrslimit.reset_index(drop=True)
corrslimit = corrslimit.loc[corrslimit['pathogen'] != 'staphylococcus_aureus', :]

corrslimit.to_csv(f'{OUT_DIR}/subset_pearson_correlations.csv', index=False)

salmdat = resdat.loc[resdat['pathogen'].isin(['salmonella_typhi', 'salmonella_paratyphi']), :]

salmdat['all_abx'] = salmdat[['pathogen', 'sample_id', 'abx_class']].groupby(
    ['pathogen', 'sample_id'])['abx_class'].transform(lambda x: ','.join(x))
reqabx = set(['aminopenicillin', 'chloramphenicol', 'sulfa', 'fluoroquinolone'])
salmdat2 = salmdat[salmdat.all_abx.str.split(',').map(reqabx.issubset)]
amino_res = list(salmdat2.loc[(salmdat2['abx_class'] == 'aminopenicillin') & (
    salmdat2['resistance'] == 'resistant'), 'sample_id'].values)
chlor_res = list(salmdat2.loc[(salmdat2['abx_class'] == 'chloramphenicol') & (
    salmdat2['resistance'] == 'resistant'), 'sample_id'].values)
sulfa_res = list(salmdat2.loc[(salmdat2['abx_class'] == 'sulfa') & (
    salmdat2['resistance'] == 'resistant'), 'sample_id'].values)


mdr_ids = set(amino_res) & set(chlor_res) & set(sulfa_res)

salmmdr = salmdat2.loc[:, ['cases', 'pathogen', 'sample_id', 'source']].drop_duplicates()

salmmdr['abx_class'] = 'mdr'

salmmdr.loc[salmmdr['sample_id'].isin(mdr_ids), 'resistance'] = 'resistant'
salmmdr.loc[~salmmdr['sample_id'].isin(mdr_ids), 'resistance'] = 'susceptible'

salmfluoro = salmdat2.loc[salmdat2['abx_class'] == 'fluoroquinolone', [
    'cases', 'pathogen', 'sample_id', 'source', 'abx_class', 'resistance']].drop_duplicates()

salmcor = salmmdr.append(salmfluoro)

salmcor.loc[salmcor['resistance'] == 'susceptible', 'resistance'] = 0
salmcor.loc[salmcor['resistance'] == 'resistant', 'resistance'] = 1

pathogen_list = []
abx_a_list = []
abx_b_list = []
nres_a_only = []
nres_b_only = []
nres_both = []
nsuscept_both = []

for bug in salmcor['pathogen'].unique():
    drugs = list(studycombos.loc[studycombos['pathogen'] == bug, 'abx_class'])
    for abxa in drugs:
        for abxb in [drug2 for drug2 in drugs if drug2 not in abxa]:
            print(datetime.now().strftime("%H:%M:%S") + bug + '-' + abxa + '-' + abxb)
            subset = salmcor.loc[(salmcor['pathogen'] == bug) & (
                salmcor['abx_class'].isin([abxa, abxb])), :]
            pivot = subset.pivot_table(
                index=['sample_id', 'pathogen', 'cases'], columns='abx_class',
                values='resistance').reset_index()
            pivot.rename(columns={abxa: 'abx_a', abxb: 'abx_b'}, inplace=True)

            pathogen_list.append(bug)
            abx_a_list.append(abxa)
            abx_b_list.append(abxb)
            nres_a_only.append(
                sum(pivot.loc[(pivot['abx_a'] > 0) & (pivot['abx_b'] == 0), 'cases']))
            nres_b_only.append(
                sum(pivot.loc[(pivot['abx_b'] > 0) & (pivot['abx_a'] == 0), 'cases']))
            nres_both.append(sum(pivot.loc[(pivot['abx_a'] > 0) & (pivot['abx_b'] > 0), 'cases']))
            nsuscept_both.append(
                sum(pivot.loc[(pivot['abx_a'] == 0) & (pivot['abx_b'] == 0), 'cases']))

coresist = pd.DataFrame({'pathogen': pathogen_list, 'abx_a': abx_a_list, 'abx_b': abx_b_list,
                         'n_res_a_only': nres_a_only, 'n_res_b_only': nres_b_only,
                         'n_res_ab': nres_both,
                         'n_suscept_ab': nsuscept_both})

coresist['prop_ab'] = coresist['n_res_ab'] / \
    (coresist['n_res_a_only'] + coresist['n_res_b_only']
     + coresist['n_res_ab'] + coresist['n_suscept_ab'])
coresist['prop_a'] = (coresist['n_res_a_only'] + coresist['n_res_ab']) / (
    coresist['n_res_a_only'] + coresist['n_res_b_only']
    + coresist['n_res_ab'] + coresist['n_suscept_ab'])
coresist['prop_b'] = (coresist['n_res_b_only'] + coresist['n_res_ab']) / (
    coresist['n_res_a_only'] + coresist['n_res_b_only']
    + coresist['n_res_ab'] + coresist['n_suscept_ab'])

coresist['pearson'] = (coresist['prop_ab'] - coresist['prop_a'] * coresist['prop_b']) / (
    np.sqrt((coresist['prop_a'] * (1 - coresist['prop_a'])
             * coresist['prop_b'] * (1 - coresist['prop_b'])))
)

coresist.to_csv(f'{OUT_DIR}/salmonella_pearson_corrs.csv')


gonorr = resdat.loc[(resdat['pathogen'] == 'neisseria_gonorrheae') & (
    resdat['abx_class'].isin(['third_gen_ceph', 'fluoroquinolone'])), :]

# limit to sample IDs with multiple tests
gonorr = gonorr.loc[:, ['pathogen', 'abx_class', 'cases',
                        'resistance', 'sample_id']].drop_duplicates()
dupsamps = gonorr.loc[gonorr['sample_id'].duplicated(), 'sample_id'].unique()

multires = gonorr.loc[gonorr['sample_id'].isin(dupsamps), :]

multires.loc[multires['resistance'] == 'susceptible', 'resistance'] = 0
multires.loc[multires['resistance'] == 'resistant', 'resistance'] = 1

pathogen_list = []
abx_a_list = []
abx_b_list = []
nres_a_only = []
nres_b_only = []
nres_both = []
nsuscept_both = []

for bug in multires['pathogen'].unique():
    drugs = ['third_gen_ceph', 'fluoroquinolone']
    for abxa in drugs:
        for abxb in [drug2 for drug2 in drugs if drug2 not in abxa]:
            print(datetime.now().strftime("%H:%M:%S") + bug + '-' + abxa + '-' + abxb)
            subset = multires.loc[(multires['pathogen'] == bug) & (
                multires['abx_class'].isin([abxa, abxb])), :]
            pivot = subset.pivot_table(
                index=['sample_id', 'pathogen', 'cases'],
                columns='abx_class', values='resistance').reset_index()
            pivot.rename(columns={abxa: 'abx_a', abxb: 'abx_b'}, inplace=True)

            pathogen_list.append(bug)
            abx_a_list.append(abxa)
            abx_b_list.append(abxb)
            nres_a_only.append(
                sum(pivot.loc[(pivot['abx_a'] > 0) & (pivot['abx_b'] == 0), 'cases']))
            nres_b_only.append(
                sum(pivot.loc[(pivot['abx_b'] > 0) & (pivot['abx_a'] == 0), 'cases']))
            nres_both.append(sum(pivot.loc[(pivot['abx_a'] > 0) & (pivot['abx_b'] > 0), 'cases']))
            nsuscept_both.append(
                sum(pivot.loc[(pivot['abx_a'] == 0) & (pivot['abx_b'] == 0), 'cases']))

coresist = pd.DataFrame({'pathogen': pathogen_list, 'abx_a': abx_a_list, 'abx_b': abx_b_list,
                         'n_res_a_only': nres_a_only, 'n_res_b_only': nres_b_only,
                         'n_res_ab': nres_both,
                         'n_suscept_ab': nsuscept_both})

coresist['prop_ab'] = coresist['n_res_ab'] / \
    (coresist['n_res_a_only'] + coresist['n_res_b_only']
     + coresist['n_res_ab'] + coresist['n_suscept_ab'])
coresist['prop_a'] = (coresist['n_res_a_only'] + coresist['n_res_ab']) / (
    coresist['n_res_a_only'] + coresist['n_res_b_only']
    + coresist['n_res_ab'] + coresist['n_suscept_ab'])
coresist['prop_b'] = (coresist['n_res_b_only'] + coresist['n_res_ab']) / (
    coresist['n_res_a_only'] + coresist['n_res_b_only']
    + coresist['n_res_ab'] + coresist['n_suscept_ab'])

coresist['pearson'] = (coresist['prop_ab'] - coresist['prop_a'] * coresist['prop_b']) / (
    np.sqrt((coresist['prop_a'] * (1 - coresist['prop_a'])
             * coresist['prop_b'] * (1 - coresist['prop_b']))))

coresist.to_csv(f'{OUT_DIR}/gonorrhea_pearson_corrs.csv')


staphbugs = list(studycombos.loc[(studycombos['pathogen'] == 'staphylococcus_aureus') & (
    studycombos['kevin_catrin_exclude'] == 0), 'abx_class'].values)

staph = resdat.loc[(resdat['pathogen'] == 'staphylococcus_aureus')
                   & (resdat['abx_class'].isin(staphbugs)), :]

# limit to sample IDs with multiple tests
staph = staph.loc[:, ['pathogen', 'abx_class', 'cases',
                      'resistance', 'sample_id']].drop_duplicates()
dupsamps = staph.loc[staph['sample_id'].duplicated(), 'sample_id'].unique()

multires = staph.loc[staph['sample_id'].isin(dupsamps), :]

multires.loc[multires['resistance'] == 'susceptible', 'resistance'] = 0
multires.loc[multires['resistance'] == 'resistant', 'resistance'] = 1

pathogen_list = []
abx_a_list = []
abx_b_list = []
nres_a_only = []
nres_b_only = []
nres_both = []
nsuscept_both = []

for bug in multires['pathogen'].unique():
    drugs = staphbugs
    for abxa in drugs:
        for abxb in [drug2 for drug2 in drugs if drug2 not in abxa]:
            print(datetime.now().strftime("%H:%M:%S") + bug + '-' + abxa + '-' + abxb)
            subset = multires.loc[(multires['pathogen'] == bug) & (
                multires['abx_class'].isin([abxa, abxb])), :]
            pivot = subset.pivot_table(
                index=['sample_id', 'pathogen', 'cases'], columns='abx_class',
                values='resistance').reset_index()
            pivot.rename(columns={abxa: 'abx_a', abxb: 'abx_b'}, inplace=True)

            pathogen_list.append(bug)
            abx_a_list.append(abxa)
            abx_b_list.append(abxb)
            nres_a_only.append(
                sum(pivot.loc[(pivot['abx_a'] > 0) & (pivot['abx_b'] == 0), 'cases']))
            nres_b_only.append(
                sum(pivot.loc[(pivot['abx_b'] > 0) & (pivot['abx_a'] == 0), 'cases']))
            nres_both.append(sum(pivot.loc[(pivot['abx_a'] > 0) & (pivot['abx_b'] > 0), 'cases']))
            nsuscept_both.append(
                sum(pivot.loc[(pivot['abx_a'] == 0) & (pivot['abx_b'] == 0), 'cases']))

coresist = pd.DataFrame({'pathogen': pathogen_list, 'abx_a': abx_a_list, 'abx_b': abx_b_list,
                         'n_res_a_only': nres_a_only, 'n_res_b_only': nres_b_only,
                         'n_res_ab': nres_both,
                         'n_suscept_ab': nsuscept_both})

coresist['prop_ab'] = coresist['n_res_ab'] / \
    (coresist['n_res_a_only'] + coresist['n_res_b_only']
     + coresist['n_res_ab'] + coresist['n_suscept_ab'])
coresist['prop_a'] = (coresist['n_res_a_only'] + coresist['n_res_ab']) / (
    coresist['n_res_a_only'] + coresist['n_res_b_only']
    + coresist['n_res_ab'] + coresist['n_suscept_ab'])
coresist['prop_b'] = (coresist['n_res_b_only'] + coresist['n_res_ab']) / (
    coresist['n_res_a_only'] + coresist['n_res_b_only']
    + coresist['n_res_ab'] + coresist['n_suscept_ab'])

coresist['pearson'] = (coresist['prop_ab'] - coresist['prop_a'] * coresist['prop_b']) / (
    np.sqrt((coresist['prop_a'] * (1 - coresist['prop_a'])
             * coresist['prop_b'] * (1 - coresist['prop_b']))))

coresist.to_csv(f'{OUT_DIR}/staph_pearson_corrs.csv')
