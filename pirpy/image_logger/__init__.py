'''
Generates a log from images inside a folder, including some default image infos
and an estimation of the FWHM.
'''

import glob
from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from photutils import find_peaks
from os.path import basename

from ..io.io_tool import path_join
from ..math.gaussian import calc_fwhm

from ..log import log

__all__ = ['ImageLogger']

default_keywords = ('BITPIX','EQUINOX', 'OBJECT', 'OBSERVER', 'TELESCOP', 'DATE-OBS', 'INSTRUME', 'EPOCH', 'FILTER')

class ImageLogger(object):
    def __init__(self, folder):
        '''
        The args are optional fields.
        '''
        self._folder = folder
        self._imlist = []
        self._fields = set([])

        self._process_folder()

    def _open_file(self, name):
        f = fits.open(name)
        return f[0].header, f[0].data

    def _find_peaks(self, data):
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        threshold = median + (10.0 * std)
        tbl = find_peaks(data, threshold)
        return tbl[:][0], tbl[:][1]

    def _get_fields(self, header, *args):
        '''
        The args are the names of the fields to get.
        '''
        result = {}
        for i in args:
            try:
                result[i] = header[i]
            except:
                result[i] = None

        return result

    def _process_file(self, name):
        f = fits.open(name)
        header = f[0].header
        data = f[0].data

        imdata = self._get_fields(header, *default_keywords)
        imdata['SIZE'] = "%ix%i" % (data.shape[0], data.shape[1])

        x, y = self._find_peaks(data)
        imdata['FWHM'] = calc_fwhm(data, x, y, max_r=30)

        imdata['FILE'] = basename(name)

        self._imlist.append(imdata)
        self._fields = self._fields.union(set(imdata.keys()))

    def _process_folder(self):
        files = glob.glob(path_join(self._folder,'*.fits'))

        for i in files:
            self._process_file(i)

    def get_table(self):
        #TODO: Sort the columns, showing first the name
        t = Table(names=self._fields, dtype=['S32']*len(self._fields))
        for i in self._imlist:
            t.add_row(i)

        return t

    def save_log(self, fname):
        ascii.write(self.get_table(), fname, format='fixed_width')
        
