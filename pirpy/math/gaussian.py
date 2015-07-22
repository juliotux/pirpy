import numpy as np
from scipy.optimize import curve_fit

from ..mp import mult_ret

def cut_region(data, xc, yc, side):
    '''
    Cuts a box region of the data, with `side` centered on (`xc`,`yc`)
    '''
    s = data.shape
    x0, x1 = max(0, xc-side), min(s[1]-1, xc+side)
    y0, y1 = max(0, yc-side), min(s[0]-1, yc+side)

    arr = data[y0:y1+1, x0:x1+1]

    return arr, x0, y0

def gaussian(x, x0, sigma, ampitude, base):
    return (amplitude/(sigma*np.sqrt(2*np.pi))*np.exp(-(x - x0)**2/(2*sigma**2))) + base

def _1d_0center_gaussian(x, sigma, amplitude, base):
    return (amplitude/(sigma*np.sqrt(2*np.pi))*np.exp(-(x - 0)**2/(2*sigma**2))) + base

def calc_r(data, xc, yc):
    '''
    This function gets a 2D array and returns the distance of each pixel to the
    (`xc`,`yc`) point and an array with the values of the original pixel at this
    distance.
    '''
    s = data.shape

    f = np.zeros(s[0]*s[1], dtype=float)
    r = np.zeros(s[0]*s[1], dtype=float)

    for y in range(s[0]):
        for x in range(s[1]):
            r[y*s[1] + x] = np.sqrt((x-xc)**2 + (y-yc)**2)
            f[y*s[1] + x] = data[y][x]

    return r, f

def _sigma(p):
    '''
    Returns: sigma
    '''
    r, f = p
    return curve_fit(_1d_0center_gaussian, r, f)[0][0]

def calc_fwhm(data, xc, yc, max_r=None, nprocess=1):
    '''
    Calculate the median fwhm of the gaussian adjust of the `data` for a given
    (`xc`,`yc`) position or list of positions. If lists of positions are given,
    you need to specify max_r.
    '''
    if isinstance(xc, (tuple, list, np.ndarray)):
        if max_r is None:
            raise ValueError('You specified a list of positions and not specified a max_r.')
        else:
            params = [None]*len(xc)
            for i in range(len(xc)):
                d, x, y = cut_region(data, xc[i], yc[i], max_r)
                x = xc[i] - x
                y = yc[i] - y
                params[i] = calc_r(d, x, y)
    else:
        params = [calc_r(data, xc, yc)]

    sigmas = mult_ret(_sigma, params, nprocess)
    return 2.355*np.nanmedian(sigmas)
