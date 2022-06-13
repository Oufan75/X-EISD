"""
This is a sample script to run eisd with N-terminal drk SH3 data that is provided in the data directory.
"""

import numpy as np
np.random.seed(91)
import os
import pandas as pd

from eisd.utils import meta_data
from eisd.parser import read_data
from eisd.utils import make_pairs
from eisd.optimizer import main

if __name__ == '__main__':

    # parameters
    data_path = "data/"         # path to experimental data and structure pools
    structure = 'ml'      # ['trades_uf', 'mixed', 'ensemble'], trades_uf refers to the Random pool
    run_mode = 'all'
    protein = 'drk'
    opt_type = 'max'    # optimization type: 'max', 'mc', None
    beta = 0.1        # hyperparameter for 'mc' opt_type (Metropolis Monte Carlo)
    
    filenames = meta_data(data_path, protein)
    exp_data = read_data(filenames['exp'], mode='exp')
    bc_data = read_data(filenames[structure], mode=structure)

    # run_mode: all
    if run_mode == 'all':
        if opt_type == 'mc':
            abs_output = "local/%s/mc_all_%s"%(structure, protein)
        else:
            abs_output = "local/%s/rl_noejc2/max_all_%s2"%(structure, protein)

        if not os.path.exists(abs_output):
            os.makedirs(abs_output)

        main(exp_data, bc_data, epochs=100, mode=['jc', 'noe', 'pre', 'fret', 'cs'], beta=beta, opt_type=opt_type, output_dir=abs_output, verbose=True) #

    # run_mode: dual
    elif run_mode == 'dual':
        pairs = [['jc', 'pre']]
        for pair in pairs:
            if opt_type == 'mc':
                abs_output = "local/%s/L+E+/mc_%s_%s_%s2"%(structure, pair[0], pair[1], protein)
            else:
                abs_output = "local/%s/max_%s_%s_%s"%(structure, pair[0], pair[1], protein)

            if not os.path.exists(abs_output):
                os.makedirs(abs_output)
            main(exp_data, bc_data, epochs=100, mode=pair, beta=beta, opt_type=opt_type, output_dir=abs_output, verbose=True)

    # run_mode: single
    elif run_mode == 'single':
        single_modes = ['cs'] #, 'fret', 'cs', 'rdc', 'rh']
        for mode in single_modes:
            if opt_type == 'mc':
                abs_output = "local/%s/L+E+/%s_%s_%s"%(structure, str(opt_type), mode, protein)
            else:
                abs_output = "local/%s/%s_%s_%s"%(structure, str(opt_type), mode, protein)

            if not os.path.exists(abs_output):
                os.makedirs(abs_output)

            main(exp_data, bc_data, epochs=100, mode=mode, beta=beta, opt_type=opt_type, output_dir=abs_output, verbose=True)

