"""
This is a sample script to run eisd with N-terminal drk SH3 data that is provided in the data directory.
"""

import numpy as np
np.random.seed(91)
import os
import pandas as pd

from eisd.parser import read_data
from eisd.utils import make_pairs
from eisd.optimizer import XEISD

if __name__ == '__main__':

    # path to experimental data and structure pools
    # supports [cs, fret, jc, rdc, rh, pre, noe, saxs] scoring and optimization
    relative_path = 'exp_data'
    exp_data_path = {
        'pre' : os.path.join(relative_path, "drksh3_pres.txt"),
        'jc'  : os.path.join(relative_path, "drksh3_JC.txt"),
        'cs'  : os.path.join(relative_path, "drksh3_CS.txt"),
    }    
    # back_calc files should have first column as index, no header, and separated by comma
    # if your files are in a different format; please adjust accordingly in eisd/parser.py
    relative_path = 'back_calc_data'
    bc_data_path =  {
        'pre' : os.path.join(relative_path, "rl_pre.txt"),
        'jc'  : os.path.join(relative_path, "rl_jc.txt"),
        'cs'  : os.path.join(relative_path, "rl_cs.txt"),
    }
    # define back calculation uncertainties
    # refer to Lincoff et al. 2020 for details
    bc_errors = {
        'pre': 0.0001,
        'noe': 0.0001,
        'saxs': 0.006,
        'fret': 0.0074,
        'rh': 0.812,
        'rdc': 0.88,
        'cs': {'C': 1.31, 'CA': 0.97, 'CB': 1.29, 'H': 0.38, 'HA': 0.29} #reported from UCBShifts
        # J-coupling errors set by default
    }

    # other parameters
    run_mode = 'all'    # supports single, dual, multiple and all data types optimization:
                        # for specific joint data combinations: use [data_type1, data_type2, ...], ex. ['jc', 'pre']
    resnum = 59         # protein residue number for SAXS calculations
    ens_size = 100      # ensemble size
    pool_size = 400     # initial conformer number
    trials = 100        # number of optimization runs 
    opt_type = 'max'    # optimization type: 'max', 'mc'
    beta = 0.1          # hyperparameter for 'mc' opt_type (Metropolis Monte Carlo)
    verbose = False     # if print warning/debugging/output infos
    
    exp_data = read_data(exp_data_path, mode='exp')
    bc_data = read_data(bc_data_path, mode='bc', bc_errors=bc_errors)
    xeisd_optimization = XEISD(exp_data, bc_data, pool_size=pool_size, nres=resnum, verbose=verbose)

    # run_mode: all
    if run_mode == 'all':
        abs_output = 'local/%s_all/'%(opt_type)    
        xeisd_optimization.optimize(trials, mode='all', ens_size=ens_size, opt_type=opt_type, output_dir=abs_output) 

    # run_mode: dual
    elif run_mode == 'duals':
        pairs = make_pairs(list(exp_data.keys()))
        #print(pairs)
        for pair in pairs:
            abs_output = 'local/%s_%s_%s/'%(opt_type, pair[0], pair[1])
            xeisd_optimization.optimize(trials, mode=pair, ens_size=ens_size, opt_type=opt_type, output_dir=abs_output)

    # run_mode: single
    elif run_mode == 'singles':
        #single_modes = ['cs', 'pre']
        for mode in exp_data.keys():
            abs_output = 'local/%s_%s/'%(opt_type, mode)
            xeisd_optimization.optimize(trials, mode=mode, ens_size=ens_size, opt_type=opt_type, output_dir=abs_output)
    
    else:
        abs_output = 'local/'      # outputs save to
        xeisd_optimization.optimize(trials, mode=run_mode, ens_size=ens_size, beta=beta, opt_type=opt_type, output_dir=abs_output)

