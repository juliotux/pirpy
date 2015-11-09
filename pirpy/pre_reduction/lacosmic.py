'''
Handle astroscrappy functions for multiple files.
'''

from astroscrappy import detect_cosmics

from ..mp import mult_ret
from ..math.ccd_type import verify_ccdlist

def __cosmics(ccddata):
    # ccddata: ccdproc.CCDData instance
    ccddata.data = detect_cosmics(ccddata.data)[1]
    return ccddata

def lacosmic_subtract(ccdlist, **kwargs):
    '''
    Subtract the cosmic rays with the L.A. Cosmic method.

    Parameters:
        ccdlist : ~ccdproc.CCDData~ or list of ~ccdproc.CCDData~
            The list of the data to be subtracted.

    Returns:
        ccdlist : ~ccdproc.CCDData~ or list of ~ccdproc.CCDData~
            The subtracted data.
    '''
    ccdlist = verify_ccdlist(ccdlist)

    return mult_ret(__cosmics, ccdlist, kwargs.get('nprocess',1))
