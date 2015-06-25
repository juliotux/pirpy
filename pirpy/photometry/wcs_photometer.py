'''
Handles the photometry functions using sky coordinates.
'''

from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from .photobject import PhotObject
from .allowed_algorithms import allowed, todo

__all__ = ['WCSPhotometer']

class WCSPhotometer(object):
    '''
    This class handles the photometry process basic functions using the
    positions in sky coordinates.
    '''
    def __init__(self, algorithm, filter=None, *args, **kwargs):
        '''
        Parameters:
            algorithm : string
                The algorithm will determine what routine will do the
                calculation of the photometry. It can assume the following
                values:
                    - 'photutils_aperture' : using the photutils astopy filiated
                                             package with aperture photometry.
                    - 'photutils_psf' : using the photutils astopy filiated
                                        package with PSF/PRF photometer. TODO
                    - 'sextractor' : using the sextractor algorithm implemented
                                     via SEP (Source Extraction and Photometry)
                                     package (aperture).
            filter : string
                The filter of the photometry, just for organization and avoid
                mistakes in the process.
        '''

        if algorithm not in allowed:
            if algorithm in todo:
                raise ValueError('The ' + algorithm + ' method is not implemented in this version. Please try another.')
            else:
                raise ValueError('The ' + algorithm + ' method is unknown. Please choose one from the list.')

        self._algorithm = algorithm
        self._filter = filter
        self._objects = []

        if self.algorithm == 'sextractor':
            from .sep_photometry import *

    def _load_image(self, fname):
        '''
        Loads one fits file, returning the WCS, header and the data matriz.
        Internal use.
        '''
