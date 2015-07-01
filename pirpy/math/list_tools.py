'''
Tools to handle single objects and lists. It's useful when you don't know if
the user is passing a single value or a list of values to a function.
'''

import numpy as np

def to_list(var):
    '''
    Check if a variable can be treated as a iterable list. If not, transform it
    to one.
    '''
    if not isinstance(var, (tuple, list, set, np.ndarray)):
            var = [var]

    return var

def match_lengths(listoflists, length):
    '''
    Check if the lists contained in a list have the same length. If the length of
    of each one is one, transform it to the givend length.
    '''
    l = listoflists
    for i in range(len(l)):
        if len(l[i]) != length and len(l[i]) == 1:
            l[i] = l[i]*length

    return l
