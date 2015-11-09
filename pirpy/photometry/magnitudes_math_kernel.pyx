import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport exp, sqrt, fabs, log10
from cython.parallel import prange

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def mag2flux(double flux, double f0=1):
    return -2.5*log10(flux/f0)

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def flux2mag(double mag, double f0=1):
    return f0 * 10**(-0.4*mag)

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef class MonteCarloMag():
    cdef readonly np.ndarray obj_fluxes, ref_fluxes, ref_mags

    def iteration(double obj_flux, double ref_flux, np.ndarray ref_mag):
        cdef double obj_mag = flux2mag(obj_flux)
        cdef int n = len(ref_flux)
        cdef double [:] results = np.zeros(n, dtype=double)

        for i in prange(n, nogil=True):
            results[i] = obj_mag + (ref_mag[i], flux2mag(ref_flux[i]))

        return np.median(results)
