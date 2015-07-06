'''
Handles the photometry functions using sky coordinates.
'''

from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from .photobject import PhotObject, PhotColection
from .allowed_algorithms import allowed, todo
from .pos_cat import PositionCatalog
from ..math.list_tools import to_list, match_lengths

from ..log import log

__all__ = ['WCSPhotometer']

class WCSPhotometer(object):
    '''
    This class handles the photometry process basic functions using the
    positions in sky coordinates.
    '''
    def __init__(self, algorithm, filter=None, log_file=None, *args, **kwargs):
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
            log_file : string
                The filename to write the log. If not set, just screen log will
                be displayed.
        '''

        if algorithm not in allowed:
            if algorithm in todo:
                raise ValueError('The ' + algorithm + ' method is not implemented in this version. Please try another.')
            else:
                raise ValueError('The ' + algorithm + ' method is unknown. Please choose one from the list.')

        self._algorithm = algorithm
        self._filter = filter
        self._objects = PhotColection(filter)

        self.position_catalog = PositionCatalog()

        self._file_queue = set([])

        if log_file is not None:
            log.enable_log_to_file(log_file)

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

    @property
    def phot_results(self):
        return self._objects

    def _load_image(self, fname):
        '''
        Loads one fits file, returning the WCS, header and the data matriz.
        Internal use.

        Returns:
            wcs : ~astropy.wcs.WCS~
                The astropy WCS from the image.
            data : ~np.ndarray~
                The data of the image.
            jd : float
                The julian date of the image.
        '''
        f = fits.open(fname)
        wcs = WCS(f[0].header)
        data = f[0].data
        jd = float(f[0].header['JD'])
        log.debug("Image %s sucessful loaded." % fname)
        return wcs, data, jd

    def _get_radec(self, x, y, wcs):
        '''
        Returns the RA and DEC coordinates from a list of (x, y) positions and a
        wcs.
        '''
        ra, dec = wcs.wcs_pix2world(x, y, 0)
        return ra, dec

    def _get_id(self, ra, dec, add_new=True):
        '''
        Returns the id of the best match from a RA,DEC pair in the position catalog.
        '''
        ra = to_list(ra)
        dec = to_list(dec)

        ids = [None]*len(ra)
        for i in range(len(ra)):
            ids[i] = self.position_catalog.match_point(ra[i], dec[i], add_new=add_new)

        return ids

    def queue_files(self, fnames):
        '''
        Queue a list of files to the photometer queue.
        '''
        fnames = to_list(fnames)
        fnames = set(fnames)
        self._file_queue = self._file_queue.union(fnames)
        log.info("%i added to file_queue." % len(fnames))

    def aperture_photometry(self, r, r_in, r_out,
                            snr, bkg_method='median',
                            elipse=False,
                            add_uid=True,
                            *args, **kwargs):
        '''
        Process the photometry for the file queue.
        '''
        #TODO: now, don't handle the elipse photometry
        #TODO: Not handle error flags at this momment.
        if self.algorithm == 'sextractor':
            from . import sep_photometry as phot

        for i in self._file_queue:
            log.debug("Meassuring photometry from %s file." % i)
            try:
                wcs, data, jd = self._load_image(i)
                bkg, rms = phot.get_background(data, bkg_method, **kwargs)
                x, y = phot.detect_sources(data, phot.get_threshold(bkg, rms, snr), **kwargs)
                log.debug("%i detected sources." % len(x))
                flux, fluxerr, flag = phot.aperture_photometry(data, x, y,
                                                               r, r_in, r_out,
                                                               elipse=False,
                                                               **kwargs)
                ra, dec = self._get_radec(x, y, wcs)
                id = self._get_id(ra, dec, add_new=add_uid)
                self._objects.add_results(jd, id, flux, fluxerr, ra, dec)
            except:
                pass

