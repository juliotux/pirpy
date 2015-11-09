'''
Handles the photometry functions using sky coordinates.
'''

from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import sep

import numpy as np

from .find import find
from .psf_fitting import psf_photometry
from .photometry import ResultStore
from ..math.list_tools import to_list, match_lengths
from ..math.gaussian import calc_fwhm

from ..log import log

__all__ = ['WCSPhotometer']

def _extract_xy(positions):
    '''
    Extract separated lists for x and y coordinates from a list of coordinates
    pairs.
    '''
    x = y = [None]*len(positions)

    for i in range(len(positions)):
        if len(positions[i]) != 2:
            raise ValueError('All the elements of the `positions` argument must' +
                             ' be a pair of 2 coordinates.')
        x[i] = positions[i][0]
        y[i] = positions[i][1]

    return x, y

def _extract_abtheta(abtheta):
    '''
    Extract separated lists for a, b and theta elipse parameters from a list of
    (a,b,theta) elements.
    '''
    a = b = theta = [None]*len(abtheta)

    for i in range(len(abtheta)):
        if len(abtheta[i]) != 3:
            raise ValueError('All the elements from abtheta argument must be' +
                             ' a list or tuple of 3 elements.')
        a[i] = abtheta[i][0]
        b[i] = abtheta[i][1]
        theta[i] = abtheta[i][2]

    return a, b, theta

def _fix_data(data):
    '''
    Fix the data matrix to a format that can be read from SEP.
    '''
    try:
        data = np.array(data)
    except:
        raise ValueError('There is a problem with your data. Please check it.')

    if data.dtype.kind is 'i':
        return data.astype('int32')
    elif data.dtype.kind is 'f':
        return data.astype('float64')
    else:
        raise ValueError('Just numerical arrays can be used in this function.')


def get_threshold(background, rms, snr):
    '''
    Determines the threshold based in the signal to noise ratio (S/N).

    Parameters:
        background : float or ~np.ndarray~
            The background of the data.
        rms : float or ~np.ndarray~
            The root main square deviation of the background.
        snr :  float
            Signal to noise ratio limmit for the detection.

    Returns:
        thresh : float or ~np.ndarray~
    '''
    #Simplified version from photutils affiliated package.
    return background + (rms*snr)


def get_background(data, bkg_method='median', *args, **kwargs):
    '''
    Determine the background image or value from the data.

    Parameters:
        data : ~numpy.ndarray~
            2D image data.
        bkg_method : string
            The method to determine the background. Can be:
                'default' : That bin the image to derive the background. Ideal
                            for non-uniform backs, but is sensible to bright
                            stars.
                'median' : use the median value of the image as background
                           estimator. Ideal for uniform backgrounds.
        **kwargs will be passed integrally to the sep functions.

    Retuns:
        bkg : float or ~numpy.ndarray~
            If the method is the 'median', it returns a float value with the
            bkg. If the method is 'default', it returns a ~numpy.ndarray~.
        rms : The root main square variance of the background.
    '''
    bkg = rms = None

    if bkg_method == 'default':
        try:
            back = sep.Background(data, **kwargs)
        except:
            data = _fix_data(data)
            back = sep.Background(data, **kwargs)
        bkg = back.back()
        rms = back.globalrms
    elif bkg_method == 'median':
        mean, bkg, rms = sigma_clipped_stats(data)
    else:
        raise ValueError('Unrecognized ' + str(bkg_method) + ' method.')

    return bkg, rms

def detect_sources(data, threshold, elipse = False, *args, **kwargs):
    '''
    Detect the sources of the image.

    Parameters:
        data : ~numpy.ndarray~
            2D image data.
        threshold : float or ~numpy.ndarray~
            The threshold of the detection.
        elipse : bool
            Tell the program if you want the elipses parameters for each object.
        **kwargs will be passed integrally to the sep functions.

    Returns:
        x, y : ~numpy.ndarray~
            The positions of the detected sources.
        a, b, theta : ~numpy.ndarray~
            The parameters of the detected elipses.
    '''
    data = _fix_data(data)
    objs = sep.extract(data, threshold, **kwargs)

    if elipse:
        return objs['x'], objs['y'], objs['a'], objs['b'], objs['theta']
    else:
        return objs['x'], objs['y']

def aperture_photometry(data, x, y, r, r_in, r_out,
                        elipse = False, abtheta = None,
                        rms = None,
                        *args, **kwargs):
    '''
    Do the aperture photometry with local sky subtraction.

    Parameters:
        data : ~numpy.ndarray~
            2D image data.
        x, y : list of list of float
            The list containing the pairs (x,y) of the objects.
        r : float
            The radius of the circular aperture to do the sum.
        r_in : float
            The internal radius of the sky annulus.
        r_out : float
            The external radius of the sky annulus.
        elipse : bool
            Tell the program if you want to do the photometry with eliptical
            aperture. If True, you have to pass the 'abtheta' argument, giving
            the elipse properties.
        **The kwargs will be passed integrally to the sep functions.

    Returns:
        flux, fluxerr, flags : ~numpy.ndarray~
            The sum of the aperture, with annulus sky subtraction, and its error.
    '''
    data = _fix_data(data)

    if elipse:
        if abtheta is None:
            raise ValueError("You must give the 'abtheta' argument if you want elipse photometry.")
        a, b, theta = _extract_abtheta(abtheta)
        return sep.sum_ellipse(data, x, y, a, b, theta, r, err=rms, bkgann=(r_in, r_out), **kwargs)
    else:
        return sep.sum_circle(data, x, y, r, err=rms, bkgann=(r_in, r_out), **kwargs)

class WCSPhotometer(object):
    '''
    This class handles the photometry process basic functions using the
    positions in sky coordinates.
    '''
    def __init__(self, filter=None, log_file=None, *args, **kwargs):
        '''
        Parameters:
            filter : string
                The filter of the photometry, just for organization and avoid
                mistakes in the process.
            log_file : string
                The filename to write the log. If not set, just screen log will
                be displayed.
        '''
        self._filter = filter
        self.results = ResultStore(filter)
        self.objects = self.results.object_catalog

        self._file_queue = set([])
        self._image_groups = {}

        if log_file is not None:
            log.enable_log_to_file(log_file)

    @property
    def filter(self):
        return self._filter

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
        x, y = wcs.wcs_world2pix(ra, dec, 0)
        return np.array(x), np.array(y)

    def _get_id(self, ra, dec, add_new=True, match_limit='1arcsec'):
        '''
        Returns the id of the best match from a RA,DEC pair in the position catalog.
        '''
        ra = to_list(ra)
        dec = to_list(dec)

        '''
        ids = np.zeros(len(ra), dtype='a32')
        for i in range(len(ra)):
            ids[i] = self.results.object_catalog.match_points(ra[i], dec[i], add_new=add_new, r_lim=match_limit)

        return ids
        '''
        return self.results.object_catalog.match_point(ra, dec, add_new=add_new, r_lim=match_limit)

    def queue_files(self, fnames, group=None):
        '''
        Queue a list of files to the photometer queue.
        '''
        fnames = set(to_list(fnames))
        self._file_queue = self._file_queue.union(fnames)
        if group is not None:
            if group in self._image_groups.keys():
                self._image_groups[group] = self._image_groups[group].union(fnames)
            else:
                self._image_groups[group] = fnames
        log.info("%i added to file_queue." % len(fnames))

    def clear_queue(self):
        '''
        Cleans the file queue.
        '''
        self._file_queue.clear()

    def aperture_photometry(self, r=None, r_in=None, r_out=None,
                            snr=20, bkg_method='median', max_r=60,
                            elipse=False, add_uid=True, objects = None,
                            xy_limits = None, nprocess=1,
                            *args, **kwargs):
        '''
        Process the photometry for the file queue.

        objects are a PhotColection with the objects to do the photometry.
        '''
        #TODO: now, don't handle the elipse photometry
        #TODO: Not handle error flags at this momment.
        if objects is not None:
            id1, ra1, dec1 = objects._id_ra_dec()

        for i in self._file_queue:
            try:
                log.info("Meassuring photometry from %s file." % i)
                wcs, data, jd = self._load_image(i)
                bkg, rms = get_background(data, bkg_method, **kwargs)

                if objects is None:
                    x, y = detect_sources(data, get_threshold(bkg, rms, snr), **kwargs)
                    x = np.array(x)
                    y = np.array(y)
                    if r is None or (r_in is None and r_out is None):
                        fwhm = calc_fwhm(data, x, y, max_r=max_r, nprocess=nprocess)
                        log.info("The median FWHM of the image is %f" % fwhm)
                    ra1, dec1 = self._get_radec(x, y, wcs)
                    id1 = self._get_id(ra1, dec1, add_new=add_uid)
                else:
                    x, y = self._get_xy(ra1, dec1, wcs)
                    if r is None or r_in is None or r_out is None:
                        x2, y2 = detect_sources(data, get_threshold(bkg, rms, snr), **kwargs)
                        fwhm = calc_fwhm(data, x2, y2, max_r=max_r, nprocess=nprocess)
                        log.info("The median FWHM of the image is %f" % fwhm)

                if xy_limits is None:
                    xy_limits = [0, data.shape[0], 0, data.shape[1]]
                t = self._check_inside_shape(x, y, xy_limits, 2*max(r_out,max_r))
                id = np.array(id1)[t]
                ra = np.array(ra1)[t]
                dec = np.array(dec1)[t]

                log.info("%i objects in image %s." % (len(id), i))

                r1, r1_in, r1_out = r, r_in, r_out
                if r is None:
                    r1 = 0.68*fwhm #optimal r for a gaussian profile
                if r_in is None:
                    r1_in = 1.5*fwhm
                if r_out is None:
                    r1_out = 2.0*fwhm

                flux, fluxerr, flag = aperture_photometry(data, x[t], y[t], r1, r1_in, r1_out, elipse=False, rms=rms, **kwargs)
                self.results.add_results(i, jd,  id, x[t], y[t], flux, fluxerr)

            except:
                log.error("Image %s failed." % i)

    def psf_photometry(self, fit_mode='radial', psf_model='gaussian',
                       snr=20, bkg_method='median', box_size=40,
                       roundlim = [-2.0, 2.0], sharplim = [-5.0, 5.0],
                       find_convolution_fwhm = 6,
                       add_uid=True, objects = None, xy_limits = None, nprocess=1,
                       mask=None, match_limit='1arcsec', *args, **kwargs):
        '''
        Process the PSF photometry in the file queue.
        '''        
        for i in self._file_queue:
            try:
                log.info("Meassuring photometry from %s file." % i)
                wcs, data, jd = self._load_image(i)
                bkg, rms = get_background(data, bkg_method)

                if objects is None:
                    sources = find(data, get_threshold(bkg, rms, snr),
                                   find_convolution_fwhm, roundlim, sharplim)
                    x = sources['x']
                    y = sources['y']

                    ra1, dec1 = self._get_radec(x, y, wcs)
                    id1 = self._get_id(ra1, dec1, add_new=add_uid, match_limit=match_limit)
                else:
                    id1, ra1, dec1 = objects._id_ra_dec()
                    x, y = self._get_xy(ra1, dec1, wcs)

                if xy_limits is None:
                    xy_limits = [0, data.shape[0], 0, data.shape[1]]
                t = self._check_inside_shape(x, y, xy_limits, box_size)
                id = np.array(id1)[t]
                ra = np.array(ra1)[t]
                dec = np.array(dec1)[t]

                flux, flux_error = psf_photometry(data, zip(x[t], y[t]), box_size=box_size,
                                                  fit_mode=fit_mode, psf_model=psf_model,
                                                  compute_errors=True)

                groups = [None]
                for j in self._image_groups.keys():
                    if i in self._image_groups[j]:
                        groups.append(j)
                self.results.add_results(i, jd, id, flux['x'], flux['y'], flux['flux'], flux_error['flux'],
                                         params=flux['parameters'], params_errors=flux_error['parameters'], groups=groups)
            except:
                log.error("Image %s failed." % i)
                #raise
