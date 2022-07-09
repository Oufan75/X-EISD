"""
This module contains miscellaneous utility functions to run the main script.
"""

def make_pairs(all):
    
    pairs=[]
    for i in range(len(all)):
        for j in range(i+1, len(all)):
            pairs.append([all[i],all[j]])
    return pairs


def modes(mode, all):
    flags = {}
    for prop in all:
        flags[prop] = False

    if mode is 'all':
        return {flag:True for flag in flags}

    elif type(mode) is list:
        for flag in mode:
            flags[flag] = True
        return flags

    elif type(mode) is str:
        flags[mode] = True
        return flags

    else:
        raise ValueError("The mode in the main function is not recognized.")
