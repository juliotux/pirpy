'''
Classify files by type and fix/add header fields.
'''

import glob
from os import path, makedirs
from astropy.io import fits

from .imdata import fits_dtype

#TODO: transfer this configure fields to a config file.
class conf(object):
    bias_keywords = set(['bias','BIAS','Bias','Zero','ZERO','zero'])
    flat_keywords = set(['flat','FLAT','Flat','FlatFiled','FLATFIELD','flatfield',
                         'FF','ff','Ff','flat_field','Flat_Field','FLAT_FIELD'])
    dark_keywords = set(['dark','DARK','Dark','Escuro','ESCURO','escuro'])
    lamp_keywords = set(['lamp'])
    obstypes = {'bias':'BIAS', 'flat':'FLAT', 'dark':'DARK', 'object':'OBJECT'}
    obstype_keyword = 'OBSTYPE'
    work_dir = '~/opd_reduction/'

def identify_file(k):
    for i in conf.bias_keywords:
        if i in k:
            return 'bias'
    for i in conf.flat_keywords:
        if i in k:
            return 'flat'
    for i in conf.dark_keywords:
        if i in k:
            return 'dark'
    for i in conf.lamp_keywords:
        if i in k:
            return 'lamp'
    return 'object'

def organize_files(folder, workdir=conf.work_dir, identifier=None):
    '''
    Organize the files inside a folder.
    '''
    if identifier is None:
        outdir = path.join(path.expanduser(workdir),'raw')
    else:
        outdir = path.join(path.expanduser(workdir),path.join(identifier,'raw'))
    try:
        makedirs(outdir)
    except:
        pass

    for i in glob.glob(path.join(folder,'*.fits')):
        try:
            f = fits.getval(i, 'OBSTYPE')
        except KeyError:
            try:
                bitpix = fits.getval(i, 'BITPIX')
                f = fits.open(i)
                k = f[0].header['OBJECT']
                t = identify_file(k)
                f[0].header[conf.obstype_keyword] = conf.obstypes[t]
                f[0].header['BITPIX'] = bitpix
                f[0].data = f[0].data.astype(fits_dtype(bitpix))
                f.writeto(path.join(outdir, path.basename(i)))
                f.close()
            except:
                print('File %s failed' % i)