'''
Wrappers the basic photometry functions from SEP (Source Extraction and Photometry)
to easy use with the photometry package.

More informations about the package: https://github.com/kbarbary/sep
'''

from astropy.stats import sigma_clipped_stats
import numpy as np

__all__ = ['get_threshold', 'get_beckground', 'detect_sources', 'aperture_photometry']

try:
    import sep
except:
    raise ImportError('The SEP package is not installed. Please install it.')

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
        return data.byteswap(True).newbyteorder()
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

    if method == 'default':
        try:
            back = sep.Background(data, **kwargs)
        except:
            data = _fix_data(data)
            back = sep.Background(data, **kwargs)
        bkg = back.back()
        rms = back.globalrms
    elif method == 'median':
        mean, bkg, rms = sigma_clipped_stats()
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
    objs = sep.extract(data, threshold, **kwargs)

    if elipse:
        return objs['x'], objs['y'], objs['a'], objs['b'], objs['theta']
    else:
        return objs['x'], objs['y']

def aperture_photometry(data, positions, r, r_in, r_out,
                        elipse = False, abtheta = None,
                        *args, **kwargs):
    '''
    Do the aperture photometry with local sky subtraction.

    Parameters:
        data : ~numpy.ndarray~
            2D image data.
        positions : list of list of float
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
    x, y = _extract_xy(positions)

    if elipse:
        if abtheta is None:
            raise ValueError("You must give the 'abtheta' argument if you want elipse photometry.")
        a, b, theta = _extract_abtheta(abtheta)
        return sep.sum_ellipse(data, x, y, a, b, theta, r, bkgann=(r_in, r_out), **kwargs)
    else:
        return sep.sum_circle(data, x, y, r, bkgann=(r_in, r_out), **kwargs)

