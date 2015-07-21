'''
Handles the photometry functions using sky coordinates.
'''

from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

import numpy as np

from .photobject import PhotObject, PhotColection
from .allowed_algorithms import allowed, todo
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
        return np.array(ra), np.array(dec)

    def _check_inside_shape(self, x, y, region, limit_radius):
        '''
        Check if the x, y coordinates are inside the values x0-x1 and y0-y1.

        region is in format [x0, x1, y0, y1]
        '''
        x0, x1, y0, y1 = region

        x = np.array(x)
        y = np.array(y)

        tx0 = x > (x0 + limit_radius)
        tx1 = x < (x1 - limit_radius)
        ty0 = y > (y0 + limit_radius)
        ty1 = y < (y1 - limit_radius)

        return tx0 & tx1 & ty0 & ty1

    def _get_xy(self, ra, dec, wcs):
        '''
        Returns the x, y coordinates from a list of ra, dec positions and a
        wcs.
        '''
        x, y = wcs.wcs_world2pix(ra, dec, 1)
        return np.array(x), np.array(y)

    def _get_id(self, ra, dec, add_new=True):
        '''
        Returns the id of the best match from a RA,DEC pair in the position catalog.
        '''
        ra = to_list(ra)
        dec = to_list(dec)

        ids = [None]*len(ra)
        for i in range(len(ra)):
            ids[i] = self._objects.match_point(ra[i], dec[i], add_new=add_new)

        return np.array(ids)

    def queue_files(self, fnames):
        '''
        Queue a list of files to the photometer queue.
        '''
        fnames = to_list(fnames)
        fnames = set(fnames)
        self._file_queue = self._file_queue.union(fnames)
        log.info("%i added to file_queue." % len(fnames))

    def clear_queue(self):
        '''
        Cleans the file queue.
        '''
        self._file_queue.clear()

    def aperture_photometry(self, r, r_in, r_out,
                            snr=20, bkg_method='median',
                            elipse=False,
                            add_uid=True,
                            objects = None,
                            xy_limits = None,
                            *args, **kwargs):
        '''
        Process the photometry for the file queue.

        objects are a PhotColection with the objects to do the photometry.
        '''
        #TODO: now, don't handle the elipse photometry
        #TODO: Not handle error flags at this momment.
        if self.algorithm == 'sextractor':
            from . import sep_photometry as phot

        if objects is not None:
            id1, ra1, dec1 = objects._id_ra_dec()

        for i in self._file_queue:
            try:
            #if True:
                log.info("Meassuring photometry from %s file." % i)
                wcs, data, jd = self._load_image(i)
                bkg, rms = phot.get_background(data, bkg_method, **kwargs)
                if objects is None:
                    x, y = phot.detect_sources(data, phot.get_threshold(bkg, rms, snr), **kwargs)
                    x = np.array(x)
                    y = np.array(y)
                    log.debug("%i detected sources." % len(x))
                    ra1, dec1 = self._get_radec(x, y, wcs)
                    id1 = self._get_id(ra1, dec1, add_new=add_uid)
                else:
                    x, y = self._get_xy(ra1, dec1, wcs)
                if xy_limits is None:
                    xy_limits = [0, data.shape[0], 0, data.shape[1]]
                t = self._check_inside_shape(x, y, xy_limits, 2*r_out)
                id = np.array(id1)[t]
                ra = np.array(ra1)[t]
                dec = np.array(dec1)[t]
                log.info("%i objects in image %s" % (len(id), i))
                flux, fluxerr, flag = phot.aperture_photometry(data, x[t], y[t],
                                                               r, r_in, r_out,
                                                               elipse=False,
                                                               **kwargs)
                self._objects.add_results(jd, id, flux, fluxerr, ra, dec)
            except:
                log.error("Image %s failed." % i)
