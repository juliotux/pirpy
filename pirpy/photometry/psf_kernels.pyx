import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport exp, sqrt, fabs, cos, sin
from cython.parallel import prange

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef class moffat:
    cdef readonly np.ndarray r, result, x, y
    cdef readonly double theta, amplitude, gamma, alpha, sky, x0, y0

    cpdef radial(self, np.ndarray[double, ndim=1] r,
                 double gamma, double alpha,
                 double amplitude, double sky):
        cdef int N = len(r)
        cdef double[:] ri = r
        cdef double[:] result = np.zeros_like(r)

        cdef int i
        for i in prange(N, nogil=True):
            result[i] = sky + amplitude*(1+(ri[i]/gamma)**2)**(-alpha)
        return result

    cpdef integrate(self, double r, double theta,
                           double amplitude, double gamma, double alpha, double sky):
        return r*amplitude*(1+(r/gamma)**2)**(-alpha)

    cpdef spatial(self, np.ndarray[double, ndim=2] x, np.ndarray[double, ndim=2] y,
                  double x0, double y0,
                  double gamma, double alpha,
                  double amplitude):
        #TODO: Parallelize the code

        '''
        cdef int N = x.shape[0]
        cdef int M = x.shape[1]

        x = x.ravel()
        y = y.ravel()

        cdef double[:] xi = x
        cdef double[:] yi = y
        cdef double[:] result = np.zeros_like(xi)

        for i in prange(N*M, nogil=True):
            result[i] = amplitude*(1 + ((xi[i] - x0)**2 + (yi[i] - y0)**2)/gamma**2)**(-alpha)
        return result
        '''

        return amplitude*(1 + ((x - x0)**2 + (y - y0)**2)/gamma**2)**(-alpha)


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef class gaussian:
    cdef readonly np.ndarray r, result, x, y
    cdef readonly double amplitude, sigma, sky, x0, y0

    cpdef radial(self, np.ndarray[double, ndim=1] r,
                 double sigma, double amplitude, double sky):
        cdef int N = len(r)
        cdef double [:] ri = r
        cdef double [:] result = np.zeros_like(r)

        cdef int i
        for i in prange(N, nogil=True):
            result[i] = sky + amplitude*exp(-0.5*(ri[i]/sigma)**2)
        return result

    cpdef spatial(self, np.ndarray[double, ndim=1] x, np.ndarray[double, ndim=1] y,
                  double x0, double y0,
                  double sigma_x, sigma_y, double theta,
                  double amplitude):

        cdef double cost2 = cos(theta)**2
        cdef double sint2 = sin(theta)**2
        cdef double sin2t = sin(2*theta)
        cdef double sigx2 = 2*sigma_x**2
        cdef double sigy2 = 2*sigma_y**2

        cdef double a = (cost2/sigx2) + (sint2/sigy2)
        cdef double b = -(sin2t/(2*sigx2)) + (sin2t/(2*sigy2))
        cdef double c = (sint2/sigx2) + (cost2/sigy2)

        cdef int N = len(x)

        cdef double [:] result = np.zeros_like(x)
        cdef double [:] xi = x - x0
        cdef double [:] yi = y - y0

        for i in range(N):
            result[i] = amplitude * exp(-(a*xi[i]**2 + 2*b*xi[i]*yi[i] + c*yi[i]**2))

        return result
        #return amplitude * np.exp(-(a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) + c*((y-y0)**2)))

