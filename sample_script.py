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
    # supports cs, fret, jc, rdc, rh, pre, noe, saxs scoring
    relative_path = 'exp_data'
    exp_data_path = {
        'pre' : os.path.join(relative_path, "drksh3_pres.txt"),
        'jc'  : os.path.join(relative_path, "drksh3_JC.txt"),
        'cs'  : os.path.join(relative_path, "drksh3_CS.txt"),
    }    
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
    ens_size = 100       # ensemble size
    pool_size = 400     # initial conformer number
    opt_type = 'max'    # optimization type: 'max', 'mc'
    beta = 0.1          # hyperparameter for 'mc' opt_type (Metropolis Monte Carlo)
    abs_output = 'local/'      # outputs save to
    if not os.path.exists(abs_output):
        os.makedirs(abs_output)
    
    exp_data = read_data(exp_data_path, mode='exp')
    bc_data = read_data(bc_data_path, mode='bc', bc_errors=bc_errors)
    xeisd_optimization = XEISD(exp_data, bc_data, pool_size=pool_size, nres=resnum)

    # run_mode: all
    if run_mode == 'all':
        xeisd_optimization.optimize(10, mode='all', ens_size=ens_size, opt_type=opt_type, output_dir=abs_output) 

    # run_mode: dual
    elif run_mode == 'duals':
        pairs = make_pairs(list(exp_data.keys()))
        #print(pairs)
        for pair in pairs:
            xeisd_optimization.optimize(100, mode=pair, ens_size=ens_size, opt_type=opt_type, output_dir=abs_output)

    # run_mode: single
    elif run_mode == 'singles':
        single_modes = ['cs', 'pre']
        for mode in single_modes:
            xeisd_optimization.optimize(100, mode=mode, ens_size=ens_size, opt_type=opt_type, output_dir=abs_output)
    
    else:
        xeisd_optimization.optimize(100, mode=run_mode, ens_size=ens_size, beta=beta, opt_type=opt_type, output_dir=abs_output)

