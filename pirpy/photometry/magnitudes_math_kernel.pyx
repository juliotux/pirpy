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
def iteration(double obj_flux, np.ndarray[double, ndim=1] ref_inst_mags, np.ndarray[double, ndim=1] ref_mag):
    cdef double obj_mag = flux2mag(obj_flux)
    cdef int n = len(inst_mags)
    cdef double [:] inst_mags = ref_inst_mags
    cdef double [:] refs = ref_mag
    cdef double [:] results = np.zeros(n, dtype=double)

    cdef int i
    for i in prange(n, nogil=True):
        results[i] = obj_mag + (refs[i] - inst_mags[i])

    return np.median(results)
