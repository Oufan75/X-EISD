import os

def meta_data(abs_path, idp):
    """
    This function provides the path to the data that is used in our papers.
    Todo: replace with a download module after publication.

    Parameters
    ----------
    abs_path: str
        The path tho the directory with the main subdirectories.
    idp: str
        The IDP name. 'drk' or 'asyn'

    Returns
    -------
    dict: dictionary of dictionaries for file names to:
        - exp: eight experimental data
        - trades: eight back calc data for Trades pool
        - tades_uf: eight back calc data for Trades pool (same as trades just wanted to embarrassingly distinguish keys)
        - mixed: eight back calc data for Mixed pool
        - ensemble: eight back calc data for Ensemble pool

    """
    # experimental data file names
    relative_path = os.path.join(abs_path, idp, "experimental_data")
    if idp == 'asyn':
        EXP_DATA_FILENAME = {
        'pre' : os.path.join(relative_path, "asyn_pres.txt"),
        'jc'  : os.path.join(relative_path, "asyn_3JHH.txt"),
        'rh'  : None,
        'rdc' : None,
        'noe' : os.path.join(relative_path, "asyn_pres.txt"),
        'fret': None,
        'cs'  : os.path.join(relative_path, "asyn_mod_CS.txt"),
        'saxs': None
    }
    else:
        EXP_DATA_FILENAME = {
        'rh'  : None,
        'rdc' : os.path.join(relative_path, "drksh3_exp_rdcs.txt"),
        'pre' : os.path.join(relative_path, "drksh3_pres.txt"),
        'noe' : os.path.join(relative_path, "8AAC_noes.txt"), #8AAC res4
        'jc'  : os.path.join(relative_path, "drksh3_JC_exp_clean.txt"), #exp_clean
        'fret': None,
        'cs'  : os.path.join(relative_path, "drksh3_CS_exp_mod.txt"),
        'saxs': os.path.join(relative_path, "unfolded_saxs_exp.txt")
    }

    # back calculated data file names
    if idp == 'asyn':
        TRADES_BC_DATA_FILENAME = None
        MIXED_BC_DATA_FILENAME = None
        relative_path = os.path.join(abs_path, idp, "back_calc_data/ML")
        ML_BC_DATA_FILENAME = {
        'rh'  : None,
        'rdc' : None,
        'pre' : os.path.join(relative_path, "rl_pre.txt"),
        'noe' : os.path.join(relative_path, "rl_pre.txt"),
        'jc'  : os.path.join(relative_path, "rl_jc.txt"),
        'fret': os.path.join(relative_path, "rl_jc.txt"),
        'cs'  : os.path.join(relative_path, "rl_cs.txt"),
        'saxs': None
    }
        relative_path = os.path.join(abs_path, idp, "back_calc_data/MCSCE")
        MCSCE_BC_DATA_FILENAME = {
	'rh'  : None,
	'rdc' : None,
	'pre' : os.path.join(relative_path, "asyn_pres.txt"),
	'noe' : os.path.join(relative_path, "drksh3_NOEs_structdists.txt"),
	'jc'  : os.path.join(relative_path, "asyn_jcs.txt"),
	'fret': os.path.join(relative_path, "asyn_jcs.txt"),
	'cs'  : os.path.join(relative_path, "asyn_CS.txt"),
	'saxs': None
     }
        
        relative_path = os.path.join(abs_path, idp, "back_calc_data/ENSEMBLE_PED")
        ENSEMBLE_BC_DATA_FILENAME = {
	'pre' : os.path.join(relative_path, "asyn_ENSEMBLE_PREs.txt"),
	'jc'  : os.path.join(relative_path, "asyn_ENSEMBLE_JC.txt"),
	'fret': os.path.join(relative_path, "asyn_ENSEMBLE_FRET.txt"),
	'cs'  : os.path.join(relative_path, "asyn_ENSEMBLE_CS.txt"),
        'rh'  : None,
	'rdc' : None,
	'noe' : os.path.join(relative_path, "asyn_ENSEMBLE_FRET.txt"),
	'saxs': None
    }
    else:
        relative_path = os.path.join(abs_path, idp, "back_calc_data/Trades")
        TRADES_BC_DATA_FILENAME = {
	'rh'  : os.path.join(relative_path, "drksh3_RH_data.txt"),
	'rdc' : os.path.join(relative_path, "drksh3_RDCs_structRDCs.txt"),
	'pre' : os.path.join(relative_path, "drksh3_PREs.txt"),
	'noe' : os.path.join(relative_path, "drksh3_NOEs.txt"),
	'jc'  : os.path.join(relative_path, "drksh3_JCs.txt"),
	'fret': os.path.join(relative_path, "drksh3_FRET_structEFF.txt"),
	'cs'  : os.path.join(relative_path, "drksh3_CS_structdata.txt"),
	'saxs': os.path.join(relative_path, "drksh3_SAXS_data.txt")
    }

        relative_path = os.path.join(abs_path, idp, "back_calc_data/mixpool")
        MIXED_BC_DATA_FILENAME = {
	'rh'  : os.path.join(relative_path, "drksh3_mixpool_RH.txt"),
	'rdc' : os.path.join(relative_path, "drksh3_mixpool_RDCs.txt"),
	'pre' : os.path.join(relative_path, "drksh3_mixpool_PREs.txt"),
	'noe' : os.path.join(relative_path, "drksh3_mixpool_NOEs.txt"),
	'jc'  : os.path.join(relative_path, "drksh3_mixpool_JC.txt"),
	'fret': os.path.join(relative_path, "drksh3_mixpool_FRET.txt"),
	'cs'  : os.path.join(relative_path, "drksh3_mixpool_CS.txt"),
	'saxs': os.path.join(relative_path, "drksh3_mixpool_SAXS.txt")
    }

        relative_path = os.path.join(abs_path, idp, "back_calc_data/ENSEMBLE_PED")
        ENSEMBLE_BC_DATA_FILENAME = {
	'rh'  : os.path.join(relative_path, "drksh3_ENSEMBLE_RH.txt"),
	'rdc' : os.path.join(relative_path, "drksh3_ENSEMBLE_RDCs.txt"),
	'pre' : os.path.join(relative_path, "drksh3_ENSEMBLE_PREs.txt"),
	'noe' : os.path.join(relative_path, "drksh3_ENSEMBLE_NOEs.txt"),
	'jc'  : os.path.join(relative_path, "drksh3_ENSEMBLE_JC.txt"),
	'fret': os.path.join(relative_path, "drksh3_ENSEMBLE_FRET.txt"),
	'cs'  : os.path.join(relative_path, "drksh3_ENSEMBLE_CS.txt"),
	'saxs': os.path.join(relative_path, "drksh3_ENSEMBLE_SAXS.txt")
    }
        relative_path = os.path.join(abs_path, idp, "back_calc_data/MCSCE")
        MCSCE_BC_DATA_FILENAME = {
	'rh'  : os.path.join(relative_path, "drksh3_RH_data.txt"),
	'rdc' : os.path.join(relative_path, "drksh3_RDCs.txt"),
	'pre' : os.path.join(relative_path, "drksh3_pres.txt"),
	'noe' : os.path.join(relative_path, "drksh3_noes.txt"),
	'jc'  : os.path.join(relative_path, "drksh3_jcs.txt"),
	'fret': os.path.join(relative_path, "drksh3_FRET.txt"),
	'cs'  : os.path.join(relative_path, "drksh3_CSs.txt"),
	'saxs': os.path.join(relative_path, "drksh3_SAXS_data.txt")
    }
        
        relative_path = os.path.join(abs_path, idp, "back_calc_data/ML")
        ML_BC_DATA_FILENAME = {
        'rh'  : None,
        'rdc' : None,
        'pre' : os.path.join(relative_path, "rl_pre.txt"),
        'noe' : os.path.join(relative_path, "rl_noe.txt"),
        'jc'  : os.path.join(relative_path, "rl_jc.txt"),
        'fret': os.path.join(relative_path, "rl_fret.txt"),
        'cs'  : os.path.join(relative_path, "rl_cs.txt"),
        'saxs': None
    }

    return {'exp':EXP_DATA_FILENAME, 'trades': TRADES_BC_DATA_FILENAME, 'trades_uf': TRADES_BC_DATA_FILENAME,
            'mixed': MIXED_BC_DATA_FILENAME, 'ensemble': ENSEMBLE_BC_DATA_FILENAME, 
            'mcsce': MCSCE_BC_DATA_FILENAME, 'ml': ML_BC_DATA_FILENAME}
