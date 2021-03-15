import numpy as np
import pandas as pd


from eisd.utils import meta_data
from eisd.parser import read_data
from eisd.scorers import *

def ensemble_scorer(property, bc, data_path='data', exp_indices=None):
    """
    This is a function to calculate EISD score for a given ensemble of conformers.
    Designed based on the current esperimental data for drkSH3 only.
    One property at a time.

    Parameters
    ----------
    property: str
        name of the experimental property:
        cs, jc, noe, fret, rh, saxs, pre, rdc

    bc: array
        2D numpy array of back-calculated data for the ensemble.
        The ensemble size is the length of array.

    data_path: str
        path to the data directory in the root directory of this repository.

    exp_indices: list, default = None
        list of indices of experimental data to include.

    Returns
    -------
    float: eisd score

    """
    filenames = meta_data(data_path)
    exp_data = read_data(filenames['exp'], mode='exp')
    bc_data = read_data(filenames['ensemble'], mode='ensemble')

    bc_data[property].data = pd.DataFrame(bc)

    if exp_indices is not None:
        exp_df = exp_data[property].data
        exp_df = exp_df.iloc[exp_indices, :]
        exp_data[property].data = exp_df

    indices = np.arange(len(bc))

    if property == 'cs':
        sse, score, _ = cs_optimization_ensemble(bc_data, exp_data, indices)

    return score

