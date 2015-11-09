'''
Classes to handle data image.
'''

from astropy.io import fits
import numpy as np

__all__ = ['ImData','fits_dtype']

nbitpix_refs = {  8 : np.uint8, #(note it is UNsigned integer)
                 16 : np.int16,
                 32 : np.int32,
                -32 : np.float32,
                -64 : np.float64}

def fits_dtype(nbitpix):
    return nbitpix_refs[nbitpix]

class ImData(object):
    def __init__(self, filename=None, header=None, data=None):
        if filename is None and (header is None or data is None):
            raise ValueError('You must specify a filename or header and data variables.')

        if filename is not None:
            bitpix = fits.getval(filename,'BITPIX')
            f = fits.open(filename)
            if len(f) > 1:
                raise ValueError('Just single HDU fits image is supported.')
            self.header = f[0].header
            self.header['BITPIX'] = bitpix
            data = f[0].data
            f.close()
            self.data = data.astype(fits_dtype(bitpix))
            del data

        else:
            self.header = header
            self.data = np.array(data)

    @property
    def shape(self):
        return self.data.shape

    @property
    def exposure(self):
        return self.header['EXPTIME']

    def set_data(self, data):
        self.data = np.array(data)

    def write_fits(filename):
        hdu = fits.PrimaryHDU(self.data, header=self.header)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(filename)