'''
Class to handle CCDData frames.
'''

from astropy.io import fits
from astropy import units as u
import numpy as np

#TODO: implement trim image

__all__ = ['CCDData', 'add_keyword', 'load_fits', 'save_fits', 'verify_ccdlist', 'set_ccd_dtype']

class CCDData(object):
    def __init__(self, data, header=None, unit=None, keywords=None, mask=None):
        #if mask is None:
        #    self._data = np.ma.masked_array(data)
        #else:
        #    self._data = np.ma.masked_array(data, mask=mask)
        self._data = np.array(data)
        self._header = header
        self._unit = unit
        if isinstance(keywords, dict):
            self._keywords = keywords
        else:
            self._keywords = dict()

    @property
    def dtype(self):
        return self._data.dtype

    @property
    def data(self):
        #return self._data.filled(fill_value=np.NaN)
        return self._data

    @data.setter
    def data(self, value):
        self._data = np.array(value)

    @property
    def keywords(self):
        return self._keywords

    @property
    def header(self):
        return self._header

    @property
    def unit(self):
        return self._unit

    def __getitem__(self, x):
        return self._data[x]

    @property
    def shape(self):
        return self.data.shape

    @property
    def exposure(self):
        return self.header['EXPTIME']

def set_ccd_dtype(ccd, dtype):
    ccd._data = ccd._data.astype(dtype)

def add_keyword(ccd, keyword, value=1, overwrite=True):
    '''
    Adds a keyword to a CCDData keywords dict.
    '''
    if keyword not in ccd.keywords.keys() or overwrite:
        ccd.keywords[keyword] = value
    #Implement a log with warn "keyword %s already exists and will not be updated"

def load_fits(fname, hdu=0, unit='adu'):
    '''
    Loads a single fits file.
    '''
    f = fits.open(fname)
    return CCDData(f[hdu].data, header=f[hdu].header, unit=unit)

def save_fits(fname, ccd, dtype=None):
    '''
    Save a CCDData frame to a fits file.
    '''
    if dtype is None:
        dtype = ccd.data.dtype
    #TODO: transform keywords to header fields.
    set_ccd_dtype(ccd, dtype)
    hdu = fits.PrimaryHDU(ccd.data, header=ccd.header)
    hdulist = fits.HDUList([hdu])
    hdu.writeto(fname)

def verify_ccdlist(ccdlist):
    '''
    Verify the ccdlist variables if it is a single CCDData instance or
    if is a list of this type of data. Returns a list with the CCDData
    '''
    if isinstance(ccdlist, list):
        for i in ccdlist:
            if not isinstance(i, CCDData):
                raise ValueError('You must give a list of CCDData' +
                                 ' instances in ccdlist variable.')
        return ccdlist

    if isinstance(ccdlist, CCDData):
        return list([ccdlist])
