from __future__ import division
import numpy as np
from scipy.optimize import curve_fit
from scipy import integrate
from astropy.nddata.utils import extract_array
from astropy.stats import sigma_clipped_stats

from . import psf_kernels

_psf_modes = set(['radial']) #TODO: spatial not working
_models = set(['gaussian', 'moffat', 'lorentzian', 'prf'])

fitting_flags = {'not_fitted':'N',             #if curve_fit couldn't fit the data
                 'sigma_discrepancy':'S',      #if the fitted sigma is discrepant from the median sigma of the image
                 'alpha_discrepancy':'A',      #if the fitted sigma is discrepant from the median sigma of the image
                 'gamma_discrepancy':'G',      #if the fitted sigma is discrepant from the median sigma of the image
                 'astrometric_error':'c',      #if the fitted coordinates differs from the original from a value higher than the errors,
                                               #or if the errors are high
                 'high_photometric_error':'p', #if the computed flux errors are high
                 }

class ModelPSF(object):
    def __init__(self, fit_mode):
        self.fit_mode = fit_mode

        if self.fit_mode == 'radial':
            self.evaluate = self.radial_kernel
            self.dtype = self.radial_dtype
            self.flag_data = self.radial_flag_data
        elif self.fit_mode == 'spatial':
            self.evaluate = self.spatial_kernel
            self.dtype = self.spatial_dtype
            self.flag_data = self.spatial_flag_data
        else:
            raise ValueError('fit_mode unrecognized.')

    @property
    def radial_dtype(self):
        raise NotImplementedError('dtype not implemented.')

    @property
    def spatial_dtype(self):
        raise NotImplementedError('dtype not implemented.')

    @staticmethod
    def radial_kernel():
        raise NotImplementedError('Radial kernel not implemented.')

    @staticmethod
    def spatial_kernel():
        raise NotImplementedError('Spatial kernel not implemented.')

    @staticmethod
    def radial_flag_data(result_array):
        raise NotImplementedError('Radial flag data not implemented.')

    @staticmethod
    def spatial_flag_data(result_array):
        raise NotImplementedError('Spatial flag data not implemented.')

    def flux_compute(self):
        raise NotImplementedError('Flux computation not implemented.')

    def fit(self):
        raise NotImplementedError('Fit not implemented.')

class GaussianPSF(ModelPSF):
    g = psf_kernels.gaussian()
    @property
    def radial_dtype(self):
        return np.dtype([('x','float64'),('y','float64'),('flux','float64'),('parameters',np.dtype(zip(['amplitude','sgima','sky'],
                        ['float64']*3)))])

    def radial_kernel(self, x, amplitude, sigma, sky):
        return self.g.radial(x, sigma, amplitude, sky)

    @staticmethod
    def radial_flag_data(result_array):
        for i in range(len(result_array)):
            result_array[i]['sigma'] = np.abs(result_array[i]['sigma'])
        median_sigma = np.nanmedian(result_array['sigma'])
        std_sigma = np.nanstd(result_array['sigma'])
        for i in range(len(result_array)):
            if np.abs(result_array[i]['sigma']-median_sigma) > 3*std_sigma:
                result_array[i]['flag'] = result_array[i]['flag'] + fitting_flags['sigma_discrepancy']

        return result_array

    @staticmethod
    def spatial_flag_data(result_array):
        return result_array

    def flux_compute(self, params, errors=None):
        if len(params) == 4:
            flux = 2*np.pi*params[0]*np.abs(params[1]*params[2])
        elif len(params) == 3:
            flux = 2*np.pi*params[0]*np.abs(params[1]*params[1])
        else:
            raise ValueError('The imput params doesn\'t correspond to a gaussian fit.')

        try:
            flux_error = np.zeros(len(params[0]))
        except TypeError:
            flux_error = 0.0

        if errors is not None:
            if len(params) == 3:
                flux_error = flux*(errors[0]/params[0] + 2*errors[1]/params[1])
            elif len(params) == 4:
                flux_error = flux*(errors[0]/params[0] + errors[1]/params[1] + errors[2]/params[2])

        return flux, flux_error

    def fit(self, x, y, data, position):
        if self.fit_mode == 'radial':
            r, f = xy2r(x, y, data, position[0], position[1])
            try:
                params, p_errors = curve_fit(self.evaluate, r, f)
            except:
                params, p_errors = np.array([np.nan]*3).reshape(3), np.array([np.nan]*9).reshape((3,3))
            p_errors = np.diag(p_errors)
            try:
                flux, flux_error = self.flux_compute(params, p_erros)
            except:
                flux, flux_error = np.nan, np.nan

            result = np.array([(flux, params)], dtype=self.dtype)
            errors = np.array([(flux_error, p_errors)], dtype=self.dtype)
            return result, errors

class MoffatPSF(ModelPSF):
    m = psf_kernels.moffat()

    @property
    def radial_dtype(self):
        return np.dtype([('x','float64'),('y','float64'),('flux','float64'),('parameters',np.dtype(zip(['amplitude','gamma','alpha','sky'],
                        ['float64']*4)))])

    def radial_kernel(self, r, amplitude, gamma, alpha, sky):
        return self.m.radial(r, gamma, alpha, amplitude, sky)

    def flux_compute(self, params, errors=None):
        r_max = params[1]*np.sqrt(10**(4/params[2]) - 1) #I(r_max) = 10^(-4)I(0)
        flux, _ = integrate.nquad(self.m.radial_integrate, [[0, r_max],[0, 2*np.pi]], args=params)
        flux_error = np.nan
        return flux, flux_error

    def fit(self, x, y, data, position):
        if self.fit_mode == 'radial':
            xp, yp = position
            r, f = xy2r(x, y, data, xp, yp)
            try:
                params, p_errors = curve_fit(self.evaluate, r, f)
                p_errors = tuple([i for i in np.diag(p_errors)])
            except:
                nantuple = tuple([np.nan]*4)
                params, p_errors = nantuple, nantuple
            try:
                flux, flux_error = self.flux_compute(params, p_errors)
            except:
                flux, flux_error = np.nan, np.nan

            result = np.array([(xp, yp, flux, params)], dtype=self.dtype)
            errors = np.array([(0, 0, flux_error, p_errors)], dtype=self.dtype)
            return result, errors

# Fitting Functions
def compute_sky(z, sigma=2, mode='median'):
    #TODO: Mode can be 'plane' too, but need to be implemented.
    mean, median, rms = sigma_clipped_stats(z)

    newz = z.ravel()
    return np.nanmedian(newz[newz < np.min(z) + sigma*rms])

def xy2r(x, y, data, xc, yc):
    r = np.sqrt((x-xc)**2 + (y-yc)**2)
    return r.ravel(), data.ravel()

def extract_data(data, indices, box_size, position):
    x, y = position
    d = extract_array(data, (box_size, box_size), (y, x))
    xi = extract_array(indices[1], (box_size, box_size), (y, x))
    yi = extract_array(indices[0], (box_size, box_size), (y, x))
    return d, xi, yi

def psf_photometry(data, positions, box_size,
                   fit_mode='radial', psf_model='gaussian',
                   compute_errors=True):
    '''
    data :
        the 2D data itself
    positions :
        the list of the positions of the stars
    box_size :
        the region around the positions to be used in the fit
    fit_mode : 'radial', 'spatial'
        the way to fit. If 'radial', the data will be modified to a 1d array,
        indexed to the distance to the star position (faster for well determined
        centroid and won't be sky subtracted, the sky will be fit). If '2d', the
        code will fit a 2D model to the original data, subtracting a
        "estimated sky".
    psf_model : 'gaussian', 'moffaf', 'lorentzian', 'prf'
        The kind of psf to be used.
    '''
    if fit_mode not in _psf_modes:
        raise ValueError('fit_mode unrecognized. The fit_mode available are:' + str(_psf_modes))
    if psf_model not in _models:
        raise ValueError('psf_model unrecognized. The available models are:' + str(_models))

    if psf_model == 'gaussian':
        psf = GaussianPSF(fit_mode)
    elif psf_model == 'moffat':
        psf = MoffatPSF(fit_mode)
    elif psf_model == 'lorentzian':
        psf = LorentzianPSF(fit_mode)

    results = np.zeros(0, dtype=psf.dtype)
    if compute_errors:
        errors = np.zeros(0, dtype=psf.dtype)

    indices = np.indices(data.shape)
    if fit_mode == 'radial':
        for i in range(len(positions)):
            d, xi, yi = extract_data(data, indices, box_size, positions[i])
            result, error = psf.fit(xi, yi, d, positions[i])
            results = np.append(results, result)
            errors = np.append(errors, error)

    try:
        results = psf.flag_data(results)
    except:
        dummy = 0

    if compute_errors:
        return results,errors
    else:
        return results