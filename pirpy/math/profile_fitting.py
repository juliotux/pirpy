from __future__ import division
import numpy as np
import warnings
from scipy.optimize import curve_fit
from astropy.modeling.parameters import Parameter
from astropy.utils.exceptions import AstropyUserWarning
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling import Fittable2DModel
from astropy.nddata.utils import extract_array, add_array, subpixel_indices

# Math Kernels
def gaussian(x, x0, sigma, amplitude, sky):
    return amplitude*np.exp(-(x - x0)**2/(2*sigma**2)) + sky

def gaussian_1d_0center(x, sigma, amplitude, sky):
    return amplitude*np.exp(-x**2/(2*sigma**2)) + sky

def moffat(x, amplitude, x0, alfa, beta, sky):
    return amplitude*(1 + ((x-x0)**2/alfa**2))**(-beta) + sky

def moffat_1d_0center(x, amplitude, alfa, beta, sky):
    return amplitude*(1 + (x**2/alfa**2))**(-beta) + sky

def lorentzian(x, amplitude, x0, gamma, sky):
    return amplitude/(1 + (x-x0)**2/gamma**2) + sky

def lorentzian_1d_0center(x, amplitude, gamma, sky):
    return amplitude/(1 + x**2/gamma**2) + sky

# Fitter Models


# Fitting functions

def calc_r(data, xc, yc):
    '''
    This function gets a 2D array and returns the distance of each pixel to the
    (`xc`,`yc`) point and an array with the values of the original pixel at this
    distance.
    '''
    s = data.shape

    f = np.zeros(s[0]*s[1])
    r = np.zeros(s[0]*s[1])

    for y from 0 <= y < s[0]:
        for x from 0 <= x < s[1]:
            r[y*s[1] + x] = np.sqrt((x-xc)**2 + (y-yc)**2)
            f[y*s[1] + x] = data[y][x]

    return r, f

def cut_region(data, xc, yc, side):
    '''
    Cuts a box region of the data, with `side` centered on (`xc`,`yc`)
    '''
    s = data.shape
    x0, x1 = max(0, xc-side), min(s[1]-1, xc+side)
    y0, y1 = max(0, yc-side), min(s[0]-1, yc+side)

    arr = data[y0:y1+1, x0:x1+1]

    return arr, xc-x0, yc-y0
