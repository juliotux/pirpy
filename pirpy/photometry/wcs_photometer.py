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

        self._file_queue = set([])

        if self.algorithm == 'sextractor':
            import .sep_photometry as phot

    @property
    def algorithm(self):
        return self._algorithm

    @property
    def filter(self):
        return self._filter

    @property
    def objects(self):
        return self._objects

    @property
    def file_queue(self):
        return self._file_queue

    def _load_image(self, fname):
        '''
        Loads one fits file, returning the WCS, header and the data matriz.
        Internal use.

        Returns:
            wcs : ~astropy.wcs.WCS~
                The astropy WCS from the image
            data : ~np.ndarray~
                The data of the image.
            jd : float
                The julian date of the image.
        '''
        f = fits.open(fname)
        wcs = WCS(f[0].header)
        data = f[0].data
        jd = float(f[0].header['JD'])
        return wcs, data, jd

    def queue_files(self, fnames):
        '''
        Queue a list of files to the photometer queue.
        '''
        if not isinstance(fnames, list):
            fnames = [fnames]

        self._file_queue = self._file_queue + set(fnames)

    def aperture_photometry(self, r, r_in, r_out):
        '''
        Process the photometry for the file queue.
        '''
        for i in self._file_queue:
            wcs, data, jd = self._load_image(i)
            #TODO: Continuar aqui
        
