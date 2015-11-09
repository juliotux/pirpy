import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport exp, sqrt, fabs
from cython.parallel import prange

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef class moffat:
    cdef readonly np.ndarray r, result
    cdef readonly double theta, amplitude, gamma, alpha, sky

    cpdef radial(self, np.ndarray[double, ndim=1] r,
                 double gamma, double alpha,
                 double amplitude, double sky):
        cdef int N = len(r)
        cdef double [:] ri = r
        cdef double [:] result = np.zeros_like(r)

        cdef int i
        for i in prange(N, nogil=True):
            result[i] = sky + amplitude*(1+(ri[i]/gamma)**2)**(-alpha)
        return result

    cpdef radial_integrate(self, double r, double theta,
                           double amplitude, double gamma, double alpha, double sky):
        return r*amplitude*(1+(r/gamma)**2)**(-alpha)


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef class gaussian:
    cdef readonly np.ndarray r, result
    cdef readonly double amplitude, sigma, sky

    cpdef radial(self, np.ndarray[double, ndim=1] r,
                 double sigma, double amplitude, double sky):
        cdef int N = len(r)
        cdef double [:] ri = r
        cdef double [:] result = np.zeros_like(r)

        cdef int i
        for i in prange(N, nogil=True):
            result[i] = sky + amplitude*exp(-0.5*(ri[i]/sigma)**2)
        return result
