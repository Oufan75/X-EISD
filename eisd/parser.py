import numpy as np
import pandas as pd

class Stack():
    def __init__(self, name, data, sigma=None, mu=None):
        self.name = name
        self.data = data
        self.sigma = sigma
        self.mu = mu



def read_data(filenames, mode, bc_errors={}):
    """
    The main function to read all the back calculated files

    Parameters
    ----------
    filenames: dict
        This parameter is a dictionary of data types with their relative path to the data file.
    bc_errors: dict
        Dict of back calculation uncertainties for each data type
    mode: str
        This parameter can be one of the following:
            - 'exp': experimental data
            - 'bc': back calculation data
            - ... or other customizable entries

    Returns
    -------
    dict: A dictionary of properties with their pandas data frame

    """
    if mode == 'exp':
        data = {}
        for key in filenames:
            # property name: Stack(name, exp data, sigma, mu)
            if key == 'fret':
                data['fret'] = Stack('fret', pd.read_csv(filenames[key]), 0.02, None)
            elif key == 'rh':
                data['rh'] = Stack('rh', pd.read_csv(filenames[key]), 0.3, None)
            else:
                data[key] = Stack(key, pd.read_csv(filenames[key]), None, None)

        return data

    elif mode == 'bc':
        data = {}
        for key in filenames:
            dfile = pd.read_csv(filenames[key], header=None, index_col=0)   
            if key == 'jc':
                data['jc'] = Stack('jc', dfile,
                        {'A': np.sqrt(0.14), 'B': np.sqrt(0.03), 'C':np.sqrt(0.08)},
                        {'A': 6.51, 'B': -1.76, 'C': 1.6})
            else:
                e = bc_errors[key]
                data[key] = Stack(key, dfile, e, None)  
                
        return data
    
 