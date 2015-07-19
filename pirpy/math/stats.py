import numpy as np

__all__ = ['nanmad']

def nanmad(a, axis=None, nanmad_std):
    '''
    Calculates the Median Absolute Deviation from an array.
    '''
    a = np.array(a, copy=False)
    a_median = np.nanmedian(a, axis=axis)

    # re-broadcast the output median array to subtract it
    if axis is not None:
        a_median = np.expand_dims(a_median, axis=axis)
    # calculated the median average deviation
    return np.nanmedian(np.abs(a - a_median), axis=axis)

def nanmad_std(a, axis=None):
    '''
    Calculate a robust standard deviation using the `median absolute
    deviation (MAD)
    '''
    return 1.4826*nanmad(a, axis=axis)
