'''
General image calibration with bias/dark/flat extraction.
'''

from astropy.io import fits
from astropy.time import Time
from os import path

#from ccdproc import image_collection

from ..ccd.ccddata import *
from ..ccd import ccdmath
from ..ccd.lacosmic import lacosmic_subtract
from ..io.io_tool import mkdir_p

class Conf(object):
    def __init__(self):
        self.ccd_dtype = 'float32'
        self.ccd_unit = 'adu'
conf = Conf()

class ImReduce(object):
    '''
    Instance for the bias/flat/dark calibration of atronomical images.

    Parameters:
        output_dir : sting
            Path to the folder to put the result images.
    '''
    def __init__(self, output_dir, **kwargs):
        self._output_dir = output_dir
        self._master_bias = None
        self._master_flat = None
        self._master_dark = None

        self._dtype = None
        self._match_keys = None

    def combine_bias(self, filelist, mean_method, **kwargs):
        '''
        Open and combines a list os files, given by filelist, to generate a master bias
        data, stored in memory.

        Parameters:
            filelist : list of string
                A list of files containig the bias frames.
            mean_method : string
                The method to use to calculate the mean of the frames to
                generate the output data. It can be 'average' or 'median':
            override (optional) : bool
                Override the actual master_bias, if it already exists.
            unit (optional) : astropy.units
                The unit of the bias images.
        '''
        #TODO: unit not yet used
        if self._master_bias != None and not kwargs.get('override',True):
            raise ValueError('Master bias already exists and won\'t be overrided.')
        else:
            self._master_bias = ccdmath.combine_bias(file_list=filelist, method=mean_method)

    #TODO: Reimplement dark correction
    #TODO: Implement default mean_method with config

    def combine_flat(self, filelist, mean_method, **kwargs):
        '''
        This routine opens a list of fits files containing the flat field frames,
        calculate the master flat frame and return it.

        Parameters:
            filelist : list of strings
                A list of files containig the flat frames.
            mean_method : string
                The method to use to calculate the mean of the frames to
                generate the output data. It can be 'average' or 'median'.
            t_exp : float
                Exposure time.
            override (optional) : bool
                Override the actual master_dark, if it already exists.
            unit (optional) : astropy.units
                The unit of the dark images.
            use_bias (optional) : bool
                Subtract the bias from all non bias corrected flat frames.
        '''
        if kwargs.get('use_bias', True):
            use_bias = kwargs.pop('use_bias', True)
            if self._master_bias == None:
                raise ValueError('No master bias found. Please, run combine_bias '+
                                 'first or set the master_bias variable by yourself.')

        if self._master_flat != None and not kwargs.get('override', True):
            raise ValueError('Master flat already exists and won\'t be overrided.')

        self._master_flat = ccdmath.combine_flat(file_list=filelist,
                                                 master_bias=self._master_bias,
                                                 skip_bias=not use_bias,
                                                 method=mean_method)

    def reduce_files(self, file_list, **kwargs):
        mkdir_p(self._output_dir)
        dtype = kwargs.get('dtype', conf.ccd_dtype)
        for i in file_list:
            print('Reducing: %s' % i)
            ccd = load_fits(i, hdu=kwargs.get('hdu', 0), unit=conf.ccd_unit)
            lacosmic_subtract(ccd)
            ccdmath.subtract_bias(ccd, self._master_bias)
            ccdmath.divide_flat(ccd, self._master_flat)
            save_fits(path.join(self._output_dir,path.basename(i)), ccd)
            del ccd
