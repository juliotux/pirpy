'''
Module to handle io operations using multiprocess or batch operations.
'''

from ccdproc import CCDData
from astropy import units as u
from os.path import join

from ..mp import mult_ret

__all__ = ['path_join','ccd_read']

################################################################################

def _path_join(path):
    '''
    Wrap the os.path.join function to work easy with arrays of file names.
    '''
    # [0] : path
    # [1] : name
    return join(path[0], path[1])

def path_join(path, filelist, **kwargs):
    '''
    Combine paths with os.path.join function, using a string or a list of
    strings as filename.

    Parameters:
        path : string
            The base path name to combine.
        filelist : string or list of strings
            The list containing the filenames to read.

    Returns:
        string or list of strings
            The data of the files.
    '''
    data = []
    if isinstance(filelist, list):
        for i in filelist:
            data.append([path, i])
        return mult_ret(_path_join, data,
                        kwargs.get('nprocess',1))
    elif isinstance(filelist, basestring):
        return _path_join([path,filelist])

################################################################################

def _read(filename, unit, dtype=None):
    '''
    Wrap the CCDData.read, including a dtype variable to
    save memory.
    '''
    dat = CCDData.read(filename, unit=unit)
    if not dtype is None:
        dat.data = dat.data.astype(dtype)
    return dat

def ccd_read(filelist, **kwargs):
    '''
    Reads .fits files into ccdproc.CCDData instances.

    Parameters:
        filelist : string or list of strings
            The list containing the filenames to read.
        unit (optional) : ~astropy.unit~
            The unit of the images.

    Returns:
        list of ~ccdproc.CCDData~
            The data of the files.
    '''
    unit = kwargs.pop('unit', u.adu)
    dt = kwargs.pop('dtype',None)
    data = []
    if isinstance(filelist, list):
        for i in filelist:
            data.append(_read(i, unit, dt))
    elif isinstance(filelist, basestring):
        data.append(_read(filelist, unit, dt))

    return data
