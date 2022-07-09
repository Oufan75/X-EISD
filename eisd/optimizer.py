import numpy as np
import pandas as pd
import os
import time

from eisd.utils import modes
from eisd.scorers import *


def monte_carlo(beta, old_total_score, new_total_score):
                
    new_probability = np.exp(beta*new_total_score)
    old_probability = np.exp(beta*old_total_score)
    # to deal with runtime error caused by large score values
    if np.any(np.isinf([old_probability, new_probability])):
        print('Runtime error... reset beta value')
        beta = 500./new_probability
        new_probability = np.exp(beta*new_score)
        old_probability = np.exp(beta*old_score)
    # accept criterion
    return np.random.random_sample() < min(1, new_probability/old_probability)


class XEISD(object):
    """
    This is the API to the X-EISD scoring to calculate and/or optimize log-likelihood of a
    disordered protein ensemble.
    Parameters
    ----------
    exp_data: dict
        Experimental data files with uncertainties.
    bc_data: dict
        Back calculation files with uncertainties.
    pool_size: int
        number of candidate conformers.
    
    verbose: bool

    """
    def __init__(self, exp_data, bc_data, nres, pool_size=None, verbose=False):
        
        self.exp_data = exp_data
        self.bc_data = bc_data
        self.verbose = verbose
        self.resnum = nres
        self.pool_size = pool_size
        if pool_size is None:
            print('Pool size not provided. Uses the JC back-calc size.')
            self.pool_size = bc_data['jc'].data.shape[0]
        
        if verbose: print("\n### Pool size: %i"%pool_size)



    def calc_scores(self, dtypes, indices=None, ens_size=100):
        '''
        Parameters
        ----------
        dtypes:
            list of data types to score 
        indices: ndarray, optional (default: None)
            This is the fastest way to get the EISD score and RMSD of selected properties for a given set of indices.
            shape: (size_of_ensemble, )
        ens_size: int
            Only used when indices not specified, to randomly select subset to score.

        Return
        ----------
        Dict of RMSDs, X-EISD scores and ensemble averages for selected conformers
        '''
        if indices is None:
            indices = np.random.choice(np.arange(self.pool_size), ens_size, replace=False)

        # initiate dict to store scores
        scores = {}
            
        for prop in self.exp_data.keys():
            scores[prop] = [0, 0, 0]
            if prop == 'jc':
                scores[prop] = [0, 0, 0, [0]]
        if 'jc' in dtypes:
            scores['jc'] = list(jc_optimization_ensemble(self.exp_data, self.bc_data, indices))
        if 'saxs' in dtypes:
            scores['saxs'] = list(saxs_optimization_ensemble(self.exp_data, self.bc_data, indices, 
                                nres=self.resnum))[:3]
        if 'cs' in dtypes:
            scores['cs'] = list(cs_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if 'fret' in dtypes:
            scores['fret'] = list(fret_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if 'noe' in dtypes:
            scores['noe'] = list(noe_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if 'pre' in dtypes:
            scores['pre'] = list(pre_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if 'rdc' in dtypes:
            scores['rdc'] = list(rdc_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if 'rh' in dtypes:
            scores['rh'] = list(rh_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
     
        return scores

            

    def optimize(self, epochs, opt_type='max', ens_size=100, mode='all', beta=0.1, 
                iters=100, output_dir=None):
        """

        Parameters
        ----------        
        opt_type: str, default 'max'
            The optimization type should be 'mc' or 'max', 'mc' for Metropolis Monte Carlo, 
            and 'max' for score maximization method.
        ens_size: int
        epochs: int
            Number of optimization trials
        mode: str or list or dict
            Data types to optimize
        beta: float
            Temperature parameter for MC optimization
        iters: int
            Number of conformer exchange attempts
        output_dir: str or path
            Directory to save optimization outputs.

        Returns
        -------

        """


        # time
        t0 = time.time()

        # switch the property
        flags = modes(mode, self.exp_data.keys())

        if output_dir is None:
            if self.verbose: print('Output directory not provided. Outputs will be saved to current directory.')
            output_dir = os.getcwd()
        elif not os.path.exists(output_dir):
            os.makedirs(output_dir)

        final_results = []
        final_indices = []
        final_best_jcoups = []

        for it in range(epochs):
            # initial scores
            indices = list(np.random.choice(np.arange(self.pool_size), ens_size, replace=False))
            old_scores = self.calc_scores([key for key in flags if flags[key]], indices)    

            new_scores = {}
            for prop in flags:
                new_scores[prop] = [0, 0, 0]
                if prop == 'jc':
                    new_scores[prop] = [0, 0, 0, [0]]
            accepted = 0
            
            for iterations in range(iters):
                pop_index = np.random.randint(0, ens_size, 1)[0]
                popped_structure = indices[pop_index]
                indices.pop(pop_index)
                struct_found = False
                while not struct_found:
                    new_index = np.random.randint(0, self.pool_size, 1)[0]
                    if new_index != popped_structure and new_index not in indices:
                        indices.append(new_index)
                        struct_found = True

                for prop in flags:
                    if flags[prop]:
                        # SAXS
                        if prop == 'saxs':
                            new_scores['saxs'] = list(saxs_optimization_ensemble(self.exp_data, 
                                                self.bc_data, None, old_scores['saxs'][2], popped_structure,
                                                new_index, ens_size, self.resnum))[:3]
                        # CS
                        if prop == 'cs':
                            new_scores['cs'] = list(cs_optimization_ensemble(self.exp_data, 
                                                self.bc_data, None, old_scores['cs'][2], popped_structure,
                                                new_index, ens_size))[:3]
                        # FRET
                        if prop == 'fret':
                            new_scores['fret'] = list(fret_optimization_ensemble(self.exp_data, 
                                                self.bc_data, None, old_scores['fret'][2], popped_structure,
                                                new_index, ens_size))[:3]
                        # JC
                        if prop == 'jc':
                            new_scores['jc'] = list(jc_optimization_ensemble(self.exp_data, 
                                                self.bc_data, None, old_scores['jc'][3], popped_structure,
                                                new_index, ens_size))
                        # NOE
                        if prop == 'noe':
                            new_scores['noe'] = list(noe_optimization_ensemble(self.exp_data, 
                                                self.bc_data, None, old_scores['noe'][2], popped_structure,
                                                new_index, ens_size))[:3]
                        # PRE
                        if prop == 'pre':
                            new_scores['pre'] = list(pre_optimization_ensemble(self.exp_data, 
                                                self.bc_data, None, old_scores['pre'][2], popped_structure,
                                                new_index, ens_size))[:3]
                        # RDC
                        if prop == 'rdc':
                            new_scores['rdc'] = list(rdc_optimization_ensemble(self.exp_data, 
                                                self.bc_data, None, old_scores['rdc'][2], popped_structure,
                                                new_index, ens_size))[:3]
                        # RH
                        if prop == 'rh':
                            new_scores['rh'] = list(rh_optimization_ensemble(self.exp_data, 
                                                self.bc_data, None, old_scores['rh'][2], popped_structure,
                                                new_index, ens_size))[:3]
            
                old_total_score = np.sum([old_scores[key][1] for key in old_scores])
                new_total_score = np.sum([new_scores[key][1] for key in new_scores])

                # optimization
                if opt_type == 'max':
                    to_accept = old_total_score < new_total_score
                elif opt_type == 'mc':
                    to_accept = monte_carlo(beta, old_total_score, new_total_score)
                else:
                    print('Opt type not supported...Abort.')
                    return 
            
                if not to_accept:
                    indices.pop(-1)
                    indices.append(popped_structure)
                else:
                    for prop in flags:
                        old_scores[prop] = new_scores[prop]

                    accepted = accepted + 1         

            s = [it, accepted]
            for prop in flags:
                # calculate scores for unoptimized data types
                if not flags[prop]:
                    if prop == 'pre':
                        old_scores['pre'][:2] = pre_optimization_ensemble(self.exp_data, self.bc_data, indices)[:2]
                    if prop == 'jc':
                        old_scores['jc'][:2] = jc_optimization_ensemble(self.exp_data, self.bc_data, indices)[:2]
                    if prop == 'cs':
                        old_scores['cs'][:2] = cs_optimization_ensemble(self.exp_data, self.bc_data, indices)[:2]
                    if prop == 'fret':
                        old_scores['fret'][:2] = fret_optimization_ensemble(self.exp_data, self.bc_data, indices)[:2]
                    if prop == 'saxs':
                        old_scores['saxs'][:2] = saxs_optimization_ensemble(self.exp_data, self.bc_data, indices,
                                                    nres=self.resnum)[:2]
                # aggregate results
                s.extend(old_scores[prop][:2])

            final_results.append(s)
            final_indices.append(indices)
            final_best_jcoups.append(old_scores['jc'][2])
            if self.verbose: print("\n### iteration: %i  (elapsed time: %f seconds)"%(it+1, time.time()-t0))

        result_header = ['index', 'accepts']
        for prop in flags:
            result_header.extend([prop+'_rmsd', prop+'_score'])
        pd.DataFrame(final_results).to_csv(os.path.join(output_dir, 'results.csv'), index=False, header=result_header)
        pd.DataFrame(final_indices).to_csv(os.path.join(output_dir, 'indices.csv'), index=False, header=False)
        if self.verbose:
            output = pd.DataFrame(final_results)
            for n in range(2, len(result_header)):
                print(result_header[n], output.iloc[:, n].mean(), output.iloc[:, n].std())
        if flags['jc']:
            pd.DataFrame(final_best_jcoups).to_csv(os.path.join(output_dir, 'best_jcoups.csv'), index=False, header=False)


