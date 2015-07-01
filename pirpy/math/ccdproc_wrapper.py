'''
Handle CCD reducing operations with multiprocessing.
'''

from ccdproc.combiner import Combiner
from ccdproc import CCDData
from ccdproc import flat_correct as fc
from ccdproc import subtract_bias as sb
from ccdproc import subtract_dark as sd
from astropy import units as u
from numpy import dtype as npdtype

from ..io import ccd_read as read
from .ccd_type import verify_ccdlist
from ..mp import mult_ret

__all__ = ['correct_flat', 'subtract_bias', 'subtract_dark', 'combine']

# Just to speedup and save memory
__main_dummy_temp_ccd = None

def __flat(ccddata):
    # [0] : CCDData to correct with __main_dummy_temp_ccd
    # [1] : kwargs
    return fc(ccddata[0], __main_dummy_temp_ccd, **ccddata[1])

def correct_flat(ccdlist, flat, **kwargs):
    '''
    Correct the flat effect of a given list of ccdproc.CCDData frames.

    Parameters:
        ccdlist : list of ccdproc.CCDData
            The list of frames to be corrected.
        flat : ccdproc.CCDData
            The master flat field frame.
        nprocess (optional) : int
            Maximum number of parallel processes. Default: 1

    Return:
        list of ccdproc.CCDData
            A list with the result frames.
    '''
    ccdlist = verify_ccdlist(ccdlist)

    global __main_dummy_temp_ccd
    __main_dummy_temp_ccd = flat

    del flat #free memory

    for i in xrange(len(ccdlist)):
        ccdlist[i] = [ccdlist[i], kwargs]

    return mult_ret(__flat, ccdlist,
                    kwargs.pop('nprocess',1))

################################################################################

def __bias(ccddata):
    # [0] : CCDData to correct with __main_dummy_temp_ccd
    return sb(ccddata, __main_dummy_temp_ccd)

def subtract_bias(ccdlist, bias, **kwargs):
    '''
    Subtract the bias effect of a given list of ccdproc.CCDData.

    Parameters:
        ccdlist : list of ccdproc.CCDData
            The list of frames to be corrected.
        bias : ccdproc.CCDData
            The master bias frame.
        nprocess (optional) : int
            Maximum number of parallel processes. Default: 1

    Return:
        list of ccdproc.CCDData
            A list with the result frames.
    '''
    ccdlist = verify_ccdlist(ccdlist)

    global __main_dummy_temp_ccd
    __main_dummy_temp_ccd = bias

    del bias #free memory

    return mult_ret(__bias, ccdlist,
                    kwargs.pop('nprocess',1))

################################################################################

def __dark(ccddata):
    # [0] : CCDData to correct with __main_dummy_temp_ccd
    # [1] : kwargs
    return sd(ccddata[0], __main_dummy_temp_ccd, **ccddata[1])

def subtract_dark(ccdlist, dark, **kwargs):
    '''
    Subtract the darkeffect of a given list of ccdproc.CCDData.

    Parameters:
        ccdlist : list of ccdproc.CCDData
            The list of frames to be corrected.
        dark : ccdproc.CCDData
            The master dark frame.
        nprocess (optional) : int
            Maximum number of parallel processes. Default: 1

    Return:
        list of ccdproc.CCDData
            A list with the result frames.
    '''
    ccdlist = verify_ccdlist(ccdlist)

    global __main_dummy_temp_ccd
    __main_dummy_temp_ccd = dark

    del bias #free memory

    for i in xrange(len(ccdlist)):
        ccdlist[i] = [ccdlist[i], kwargs]

    return mult_ret(__dark, ccdlist,
                    kwargs.pop('nprocess',1))

################################################################################

def combine(input_list, mean_method):
    '''
    Combines a list of ~ccdproc.CCDData~ into a single, using ~ccdproc.Combiner~.

    Parameters:
        input_list : list of ccdproc.CCDData
            The CCDData to be combined.
        mean_method: string
            The method to calculate the mean. Can be:
                'average'
                'median'

    Returns:
        ~ccdproc.CCDData~
            The combined data.
    '''
    comb = Combiner(input_list)

    if mean_method == 'median':
        return comb.median_combine()
    elif mean_method == 'average':
        return comb.average_combine()
