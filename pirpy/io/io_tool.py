'''
Module to handle io operations using multiprocess or batch operations.
'''

from astropy import units as u
from os.path import join
from os import makedirs
import os
import errno
from astropy.io import fits

from ..mp import mult_ret
from ..ccd.ccddata import load_fits, set_ccd_dtype

__all__ = ['path_join', 'mkdir_p']

def mkdir_p(fname):
    try:
        makedirs(fname)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(fname):
            #TODO: implement a log.info here
            pass

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
