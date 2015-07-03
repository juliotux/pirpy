'''
General image calibration with bias/dark/flat extraction.
'''

from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from ccdproc import CCDData
from ccdproc.combiner import Combiner
import numpy as np
import bottleneck as bn
from os import path

from ..math import *
from ..io.mkdir_p import mkdir_p
from ..io.io_tool import ccd_read
from lacosmic import lacosmic_subtract

class ImReduce():
    '''
    Instance for the bias/flat/dark calibration of atronomical images.

    Parameters:
        output_dir : sting
            Path to the folder to put the result images.
    '''
    def __init__(self, output_dir, **kwargs):
        self.output_dir = output_dir
        self.master_bias = None
        self.master_flat = None
        self.master_dark = None

        self.dark_exptime = 0.0
        self.flat_exptime = 0.0

        self.dtype = None
        self.match_keys = None

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
            t_exp (optional) : float
                Exposure time.
            unit (optional) : astropy.units
                The unit of the bias images.
        '''
        if self.master_bias != None and not kwargs.get('override',True):
            print('Master bias already exists and won\'t be overrided.')

        bias = Combiner(ccd_read(filelist, unit=kwargs.get('unit', u.adu)))

        if mean_method == 'average':
            self.master_bias = bias.average_combine()
        elif mean_method == 'median':
            self.master_bias = bias.median_combine()

    def combine_dark(self, filelist, mean_method, exp_time, **kwargs):
        '''
        Combine a list of files to generate the master dark frame.

        Parameters:
            filelist : list of string
                A list of files containig the dark frames.
            mean_method : string
                The method to use to calculate the mean of the frames to
                generate the output data. It can be 'average' or 'median'.
            t_exp : float
                Exposure time.
            override (optional) : bool
                Override the actual master_dark, if it already exists.
            unit (optional) : astropy.units
                The unit of the dark images.
        '''
        if kwargs.get('use_bias',True):
            if self.master_bias == None:
                raise ValueError('No master bias found. Please, run combine_bias '+
                                 'first or set the master_bias variable by yourself')
        if self.master_dark != None and not kwargs.get('override',True):
            print('Master dark already exists and won\'t be overrided.')

        self.dark_exptime = t_exp

        dark = None
        if kwargs.get('use_bias',True):
            dark = Combiner(subtract_bias(
                                ccd_read(filelist,
                                    unit=kwargs.get('unit', u.adu),
                                    dtype=dtype),
                                self.master_bias,
                                nprocess=kwargs.get('nprocess',1)))

        if mean_method == 'average':
            self.master_dark = dark.average_combine()
        elif mean_method == 'median':
            self.master_dark = dark.median_combine()

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
        '''
        if kwargs.get('use_bias',True):
            if self.master_bias == None:
                raise ValueError('No master bias found. Please, run combine_bias '+
                                 'first or set the master_bias variable by yourself.')

        if kwargs.get('use_dark',False):
            if self.master_dark == None:
                raise ValueError('No master dark found. Please, run combine_dark '+
                                 'first or set the master_dark variable by yourself.')
            if not 't_exp' in kwargs.keys():
                raise ValueError('To use dark, you must set the flat field exposure.')

        if self.master_flat != None and not kwargs.get('override',True):
            print('Master flat already exists and won\'t be overrided.')

        dtype = kwargs.get('dtype',None)
        flat = None
        if kwargs.get('use_bias',True) and kwargs.get('use_dark',False):
            flat = Combiner(subtract_dark(
                                subtract_bias(
                                    ccd_read(filelist,
                                        unit=kwargs.get('unit', u.adu)),
                                    self.master_bias,
                                    nprocess=kwargs.get('nprocess',1),
                                    dtype=dtype),
                                self.master_dark,
                                dark_exposure=self.dark_exptime,
                                data_exposure=self.flat_exptime,
                                scale=True,
                                nprocess=kwargs.get('nprocess',1)))
        elif kwargs.pop('use_bias',True):
            flat = Combiner(subtract_bias(
                                ccd_read(filelist,
                                    unit=kwargs.get('unit', u.adu),
                                    dtype=dtype),
                                self.master_bias,
                                nprocess=kwargs.get('nprocess',1)))

        if mean_method == 'average':
            flat = flat.average_combine()
        elif mean_method == 'median':
            flat = flat.median_combine()

        self.master_flat = flat.divide(bn.nanmean(flat.data)) #Normalized flat
        del flat

    def reduce_files(self, filelist, **kwargs):
        mkdir_p(self.output_dir)
        dtype = kwargs.get('dtype','float32')
        arglist = []
        for i in filelist:
            print('Reducing: ' + i)
            ccd = CCDData.read(i, unit='adu')
            ccd = lacosmic_subtract(ccd)
            ccd = subtract_bias(ccd, self.master_bias)
            ccd = correct_flat(ccd, self.master_flat)
            ccd = ccd[0]
            ccd.data = ccd.data.astype(dtype)
            ccd = ccd.to_hdu()
            ccd.writeto(path.join(self.output_dir,path.basename(i)))
            del ccd
