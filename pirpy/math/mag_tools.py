import numpy as np

__all__ = ['flux2mag','mag2flux']

def flux2mag(flux):
    return -2.5*np.log10(flux)

def mag2flux(mag):
    return 10**(-0.4*mag)
