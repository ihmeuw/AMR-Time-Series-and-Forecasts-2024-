import sys
import pandas as pd
import numpy as np
import sys
import os
import argparse
import getpass
user = getpass.getuser()
simbinpath = 'FILEPATH'.format(user)
if not simbinpath in sys.path:
    sys.path.append('FILEPATH'.format(user))
amrpath = 'FILEPATH'.format(user)
if not amrpath in sys.path:
    sys.path.append('FILEPATH'.format(user))
from amr_prep.utils.amr_io import get_amr_data
from mcod_prep.utils.mcause_io import get_mcause_data

import itertools
from typing import Dict
import numpy as np
from numpy import ndarray
from scipy.optimize import LinearConstraint, minimize
import numpy as np
from numpy import ndarray

def sample_simplex(ndim: int, size: int = 1) -> ndarray:
    # special case when n == 1
    if ndim == 1:
        return np.ones((size, ndim))

    points = np.hstack([np.zeros((size, 1)),
                        np.random.rand(size, ndim - 1),
                        np.ones((size, 1))])
    points.sort(axis=1)
    samples = np.diff(points)

    return samples


def adjust_prob_mat(prob_mat: ndarray,
                    weights: ndarray = None,
                    x0: ndarray = None,
                    rtol: float = 1e-4,
                    atol: float = 1e-4,
                    options: Dict = None) -> ndarray:
    ndim = prob_mat.shape[0]
    size = 2**ndim

    weights = np.ones((ndim, ndim)) if weights is None else weights
    A = np.empty(shape=(0, size))
    b = np.empty(0)
    w = np.empty(0)

    ind = [[0, 1]]*ndim
    for i in range(ndim):
        for j in range(i, ndim):
            ind[i] = [1]
            ind[j] = [1]
            row = np.zeros((2,)*ndim)
            row[tuple(np.meshgrid(*ind))] = 1.0
            A = np.vstack([A, row.ravel()])
            b = np.hstack([b, prob_mat[i, j]])
            w = np.hstack([w, weights[i, j]])
            ind[i] = [0, 1]
            ind[j] = [0, 1]

    # create optimization problem
    def objective(x: ndarray) -> float:
        return 0.5*w.dot((A.dot(x) - b)**2)

    def gradient(x: ndarray) -> ndarray:
        return (A.T*w).dot(A.dot(x) - b)

    def hessian(x: ndarray) -> ndarray:
        return (A.T*w).dot(A)

    if x0 is None:
        x0 = np.ones(size)/size
    bounds = np.hstack([np.zeros((size, 1)), np.ones((size, 1))])
    constraints = [LinearConstraint(np.ones((1, size)), np.ones(1), np.ones(1))]

    result = minimize(
        objective, x0,
        method="trust-constr",
        jac=gradient,
        hess=hessian,
        constraints=constraints,
        bounds=bounds,
        options=options
    )

    x = result.x
    adjusted_b = A.dot(x)
    adjusted_prob_mat = prob_mat.copy()
    if not np.allclose(adjusted_b, b, rtol=rtol, atol=atol):
        k = 0
        for i in range(ndim):
            for j in range(i, ndim):
                adjusted_prob_mat[i, j] = adjusted_b[k]
                k += 1
    return adjusted_prob_mat, x

IN_DIR = "FILEPATH"

corrs = pd.read_csv(
    f'{IN_DIR}/subset_pearson_correlations.csv')
corrs_salmonella = pd.read_csv(
    f'{IN_DIR}/salmonella_pearson_corrs.csv')
corrs = corrs.append(corrs_salmonella)

studycombos = pd.read_csv('FILEPATH'.format(user), encoding = 'latin-1')

PATHOGEN = str(sys.argv[1])
LOCATION_ID = str(sys.argv[2])

# initialize blank DataFrames
resdf = pd.DataFrame()
stgprres = pd.DataFrame()
# pull out correlations by pathogen
bugdf = corrs.loc[corrs['pathogen'] == PATHOGEN, :]
# initialize matrix of across antibiotic combinations
mtx = bugdf.pivot_table(index='abx_a', columns='abx_b', values='pearson')
bugdrugs = bugdf['abx_a'].sort_values().unique()

# create frame of STGPR draws (1k) for each antibiotic at specified location
for drug in bugdrugs:
    bugfile = pd.read_csv(
        (
            "FILEPATH"
            + studycombos.loc[
                (studycombos['pathogen'] == PATHOGEN) & (studycombos['abx_class'] == drug),
                'resistance_run_id'
            ].astype('int').astype('str')
            + "/draws_temp_0/" + str(LOCATION_ID) + ".csv"
        ).values[0]
    )
    bugfile['antibiotic'] = drug
    stgprres = stgprres.append(bugfile)

# matrix creation and sampling for ea draw
for y in range(1990, 2024):
    for draw in ["draw_" + str(i) for i in range(0, 100)]:
    # fill out diagonals from STGPR results
        for drug in bugdrugs:
            mtx.loc[drug, drug] = stgprres.loc[(stgprres['antibiotic'] == drug) & (stgprres.year_id == y), draw].values[0]
    # populate off diagonals with expected A&B prevalence based on Pearson correlations
        for drug1 in bugdrugs:
            for drug2 in [drug2 for drug2 in bugdrugs if drug2 not in drug1]:
                pearson = bugdf.loc[(bugdf['abx_a'] == drug1) & (
                    bugdf['abx_b'] == drug2), 'pearson'].values
                mtx.loc[drug1, drug2] = pearson * np.sqrt(
                    mtx.loc[drug1, drug1] * (1 - mtx.loc[drug1, drug1])
                    * mtx.loc[drug2, drug2] * (1 - mtx.loc[drug2, drug2])
                ) + mtx.loc[drug1, drug1] * mtx.loc[drug2, drug2]

    # create and adjust probability matrix
        probab_mat = mtx.values
        adjusted_prob_mat, prob = adjust_prob_mat(
            prob_mat=probab_mat, x0=sample_simplex(2**probab_mat.shape[0], size=1)[0])
        prob = prob.reshape(np.repeat([2], len(bugdrugs)))


    # initialize frame of all combinatorics
        lst = list(itertools.product([0, 1], repeat=len(bugdrugs)))
        combos = pd.DataFrame.from_records(lst)
        combos.columns = bugdrugs
    # impute proportion of each combinatoric based on binomial sampling
        for i in range(0, len(combos)):
            combos.loc[i, 'combinatoric_prop'] = prob[tuple(
                combos.loc[i, bugdrugs].astype('int').values)]
    # create ID for each combinatoric
        combos.reset_index(inplace=True)
        combos['index'] = PATHOGEN + '-' + combos['index'].astype('str')
        combos = combos.rename(columns={'index': 'combinatoric_id'})
    # reshape to long
        melted = combos.melt(id_vars=['combinatoric_id', 'combinatoric_prop'],
                         var_name='abx_class', value_name='resistant')
        melted['draw'] = draw
        melted['year_id'] = y
        melted['location_id'] = LOCATION_ID
        melted['pathogen'] = PATHOGEN
        resdf = resdf.append(melted)
# reshape to wide on draw
    resdf2 = resdf.pivot_table(index=['combinatoric_id', 'pathogen', 'location_id','year_id',
                                      'abx_class', 'resistant'],
                               columns='draw', values='combinatoric_prop').reset_index()
    outdir = 'FILEPATH' + \
        str(PATHOGEN) + '/'
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    resdf2.to_csv(outdir + str(LOCATION_ID) + '.csv', index=False)
