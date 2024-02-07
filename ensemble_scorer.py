"""
This is a sample script to score an ensemble pool.
"""

import numpy as np
np.random.seed(91)
import os
import pandas as pd

from eisd.parser import read_data
from eisd.optimizer import XEISD

if __name__ == '__main__':

    # path to experimental data and structure pools
    # supports cs, fret, jc, rdc, rh, pre, noe, saxs scoring
    relative_path = 'exp_data'
    exp_data_path = {
        'noe' : os.path.join(relative_path, "ab40_noes.txt"),
        #'pre' : os.path.join(relative_path, "asyn_pres.txt"),
        'jc'  : os.path.join(relative_path, "ab40_jc.txt"),
        'fret'  : os.path.join(relative_path, "ab40_efret.txt"),
        #'saxs': os.path.join(relative_path, "ab40_saxs.txt"),
        'cs': os.path.join(relative_path, "ab40_cs.txt"),
    }    
    relative_path = '../IDPdiff/local/esm_200T/epoch_170_sample/ab40'
    bc_data_path =  {
        'noe' : os.path.join(relative_path, "noe.txt"),
        #'pre' : os.path.join(relative_path, "pre.txt"),
        'jc'  : os.path.join(relative_path, "jc.txt"),
        'fret'  : os.path.join(relative_path, "fret.txt"),
        #'saxs': os.path.join(relative_path, "saxs.txt"),
        'cs': os.path.join(relative_path, "cs.txt"),
    }
    # define back calculation uncertainties
    bc_errors = {
        'pre': 0.1,
        'noe': 0.1,
        'saxs': 0.006,
        'fret': 0.0074,
        'rh': 0.812,
        'rdc': 0.88,
        'cs': {'C': 1.31, 'CA': 0.97, 'CB': 1.29, 'H': 0.38, 'HA': 0.29, 'N': 2.16} #reported from UCBShifts
        # J-coupling errors set by default
    }

    # other parameters
    resnum = 59                  # protein residue number for SAXS calculations
    ens_size = 100               # ensemble size
    pool_size = 200              # total conformer number
    trials = 10                 # scoring attempt number
    data_types = ['jc', 'noe', 'fret', 'cs']  # list of data types to score
    per_atom_CS_error = True
    
    exp_data = read_data(exp_data_path, mode='exp')
    bc_data = read_data(bc_data_path, mode='bc', bc_errors=bc_errors)
    xeisd_optimization = XEISD(exp_data, bc_data, pool_size=pool_size, nres=resnum)

    abs_output = '../IDPdiff/local/esm_200T/epoch_58_sample/drkn'        # outputs save to
    if not os.path.exists(abs_output):
        os.makedirs(abs_output)
    # if first time calculate
    output = pd.DataFrame(np.zeros((trials, 2*len(exp_data) + 6 if per_atom_CS_error else 2*len(exp_data))))
    # else: read existing file
    #output = pd.read_csv(os.path.join(abs_output, 'scores.csv'), index_col=0)
    
    for n in range(trials):
        # if read indices from file, else randomly select from pool
        results = xeisd_optimization.calc_scores(data_types, ens_size=ens_size, expand_CS_errors=per_atom_CS_error)
        st = 0
        for prop in data_types:
            output.iloc[n, st] = results[prop][0]
            output.iloc[n, st+1] = results[prop][1]
            st += 2
        if per_atom_CS_error:
            for atom in results['cs_per_atom_rmsd']:
                output.iloc[n, st] = results['cs_per_atom_rmsd'][atom]
                st += 1


    result_header = []
    for prop in data_types:
        result_header.extend([prop+'_rmsd', prop+'_score'])
    if per_atom_CS_error:
        result_header.extend([f'cs_{a}_rmsd' for a in results['cs_per_atom_rmsd']])
    print('{0: >35} {1: >25}'.format('mean', 'stdev'))
    for n in range(len(result_header)):
        print(f'{result_header[n]: <10} {output.iloc[:, n].mean(): >25} {output.iloc[:, n].std(): >25}')
    #output.to_csv(os.path.join(abs_output, 'scores.csv'), header=result_header)

