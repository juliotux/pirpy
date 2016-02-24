'''
General image calibration with bias/dark/flat extraction.
'''

from astropy.io import fits
from astropy.time import Time
from os import path

#from ccdproc import image_collection

from ..ccd.ccddata import *
from ..ccd.ccdmath import *
from .lacosmic import lacosmic_subtract

class ImReduce():
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
        if self.master_bias != None and not kwargs.get('override',True):
            raise ValueError('Master bias already exists and won\'t be overrided.')
        else:
            self._master_bias = combine_bias(filelist, mean_method)

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
        if kwargs.get('use_bias', True):
            if self._master_bias == None:
                raise ValueError('No master bias found. Please, run combine_bias '+
                                 'first or set the master_bias variable by yourself.')

        if self._master_flat != None and not kwargs.get('override', True):
            raise ValueError('Master flat already exists and won\'t be overrided.')

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
