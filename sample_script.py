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
    # exp data should contain index, value, error as header
    relative_path = '4ebp'
    exp_data_path = {
        "fret": os.path.join(relative_path, "exp.fret"),
        "rh": os.path.join(relative_path, "exp.rh"),
        #'pre' : os.path.join(relative_path, "drksh3_exp_pres.txt"),
        #'jc'  : os.path.join(relative_path, "drksh3_exp_JC.txt"),
        #'cs' : os.path.join(relative_path, "drksh3_exp_CS.txt"),
    }    
    # back_calc files should have first column as index, no header, and separated by comma
    # each row is a conformer and each column is back calculation aligned with the exp data
    # if your files are in a different format; please adjust accordingly in eisd/parser.py
    relative_path = '4ebp'
    bc_data_path =  {
        "fret": os.path.join(relative_path, "bc.fret"),
        "rh": os.path.join(relative_path, "bc.rh"), 
        #'pre' : os.path.join(relative_path, "drksh3_pres.txt"),
        #'jc'  : os.path.join(relative_path, "drksh3_jcs.txt"),
        #'cs'  : os.path.join(relative_path, "drksh3_CSs.txt"),
    }
    # define back calculation uncertainties
    # refer to Lincoff et al. 2020 for details
    bc_errors = {
        'pre': 0.0001,
        'noe': 0.0001,
        'saxs': 0.006,
        'fret': 0.03,
        'rh': 0.812,
        'rdc': 0.88,
        'cs': {'C': 1.31, 'CA': 0.97, 'CB': 1.29, 'H': 0.38, 'HA': 0.29} #reported from UCBShifts
        # J-coupling errors set by default
    }
    
    # other parameters
    run_mode = "all"    # supports single, dual, multiple and all data types optimization:
                        # for specific joint data combinations: use [data_type1, data_type2, ...], ex. ['jc', 'pre']
    trials = 5
    resnum = 59         # residue number only used for SAXS
    ens_size = 100      # final ensemble size
    pool_size = 10000     # initial conformer number
    opt_type = 'max'    # optimization type: 'max', 'mc'
    beta = 1         # hyperparameter for 'mc' opt_type (Metropolis Monte Carlo)
    verbose = True
    abs_output = '4ebp/max_all'      # outputs save to
    if not os.path.exists(abs_output):
        os.makedirs(abs_output)
    
    exp_data = read_data(exp_data_path, mode='exp')
    bc_data = read_data(bc_data_path, mode='bc', bc_errors=bc_errors)
    xeisd_optimization = XEISD(exp_data, bc_data, pool_size=pool_size, nres=resnum, verbose=verbose)

    # run_mode: all
    if run_mode == 'all':
        xeisd_optimization.optimize(trials, mode='all', ens_size=ens_size, beta=beta, opt_type=opt_type, output_dir=abs_output) 

    # run_mode: dual
    elif run_mode == 'duals':
        pairs = make_pairs(list(exp_data.keys()))
        for pair in pairs:
            abs_output = '%s/%s_%s_%s/'%(abs_output, opt_type, pair[0], pair[1])
            xeisd_optimization.optimize(trials, mode=pair, ens_size=ens_size, beta=beta, opt_type=opt_type, output_dir=abs_output)

    # run_mode: single
    elif run_mode == 'singles':
        for mode in exp_data.keys():
            abs_output = '%s/%s_%s/'%(abs_output, opt_type, mode)
            xeisd_optimization.optimize(trials, mode=mode, ens_size=ens_size, beta=beta, opt_type=opt_type, output_dir=abs_output)
    
    else:
        xeisd_optimization.optimize(trials, mode=run_mode, ens_size=ens_size, beta=beta, opt_type=opt_type, output_dir=abs_output)

