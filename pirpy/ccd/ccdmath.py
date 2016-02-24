import numpy as np
from .ccddata import CCDData, load_fits, add_keyword
from .lacosmic import lacosmic_subtract

try:
    import bottleneck as bn
    nmedian = bn.nanmedian
    nmean = bn.nanmean
    nstd = bn.nanstd
    nmin = bn.nanmin
    nmax = bn.nanmax
    nvar = bn.nanvar
except NameError:
    nmedian = np.nanmedian
    nmean = np.nanmean
    nstd = np.nanstd
    nmin = np.nanmin
    nmax = np.nanmax
    nvar = np.nanvar

def func_list(func, ccd_list, *args, **kwargs):
    '''
    Apply a function to a list (or set) of CCDData elements.
    '''
    if isinstance(ccd_list, (tuple, list, set)):
        for i in ccd_list:
            if not isinstance(i, CCDData):
                raise ValueError('ccd_list must contain just CCDData objects.')
    else:
        raise ValueError('ccd_list must be a list of CCDData elements.')

    ndata = len(ccd_list)

    #TODO: implement header here, or keywords
    #Now we leave the dtype free. Think more about it.
    return CCDData(func(np.array([ccd_list[i].data for i in range(ndata)]), axis=0, *args, **kwargs))#.astype(ccd_array[0].dtype))

def func_ccd(func, ccd):
    '''
    Apply a function to a single CCDData instance.
    '''
    if not isinstance(i, CCDData):
        raise ValueError('ccd argument must be a CCDData instance')

    return func(ccd.data)

def median_list(ccd_list):
    return func_list(nmedian, ccd_list)

def mean_list(ccd_list):
    return func_list(nmean, ccd_list)

#TODO: to implement with masked arrays
#def sigclipavg_list(ccd_array):
#    func_list(sigclipavg, ccd_array)

valid_methods = set(['median', 'mean'])
combine_funcs = {'median': median_list,
                 'mean': mean_list}

################################################################################

def ccd_math(op, ccd, value=None):
    '''
    Apply math operations to a CCDData.
    '''
    if isinstance(value, CCDData):
        v = value.data
    else:
        v = value

    if op in set(['+', 'add', 'sum', 'plus']):
        ccd.data = ccd.data + v
    elif op in set(['-', 'subtract', 'minus']):
        ccd.data = ccd.data - v
    elif op in set(['*', 'multiply']):
        ccd.data = ccd.data * v
    elif op in set(['/', 'divide']):
        ccd.data = ccd.data / v
    elif op == 'sqrt':
        ccd.data = np.sqrt(ccd.data)
    elif op == 'square':
        ccd.data = ccd.data**2

################################################################################

def files_or_ccd(**kwargs):
    if 'file_list' in kwargs.keys() and 'ccd_list' in kwargs.keys():
        raise ValueError('Only one of this arguments is supported at same time: file_list or ccd_list')

    if 'file_list' in kwargs.keys():
        file_list = kwargs.pop('file_list')
        ccd_list = []
        if isinstance(file_list, (tuple, list, set)):
            for i in file_list:
                ccd_list.append(load_fits(i))
        else:
            raise ValueError('file_list have a wrong format.')

    elif 'ccd_list' in kwargs.keys():
        ccd_list = kwargs.pop('ccd_list')
        if isinstance(ccd_list, (tuple, list, set)):
            for i in ccd_list:
                if not isinstance(i, CCDData):
                    raise ValueError('ccd_list must contain just CCDData objects.')
        else:
            raise ValueError('ccd_list must be a list of CCDData elements.')
    else:
        raise ValueError('You must specify a ccd_list or a file_list for this function.')

    return ccd_list

def get_method(**kwargs):
    if 'method' in kwargs.keys():
        method = kwargs.pop('method')
        if method not in valid_methods:
            raise ValueError('Unrecognized mean method.')
    else:
        method = 'median'
    #TODO: if implement masked arrays and sigma_clip statistics, load sigma_clip here.
    return method

################################################################################
#TODO: implement a "check_header_keywords" method for the following defs
#TODO: implement a function to check if the images are already corrected by bias/flat/dark
#TODO: implement a function to wirte keywords for master_bias, master_flat, or data corrected by bias/flat/dark
def combine_bias(**kwargs):
    '''
    Combine bias images in a master_bias.
    '''
    ccd_list = files_or_ccd(**kwargs)
    method = get_method(**kwargs)

    if kwargs.get('lacosmic', False):
        lacosmic_subtract(ccd_list)

    func = combine_funcs[method]
    result = func(ccd_list)
    add_keyword(result, 'master_bias', 1)
    add_keyword(result, 'combine_bias_method', method)
    return result

def subtract_bias(ccd, master_bias, skip_keywords=False):
    '''
    Subtract the master bias file from a ccd image.
    '''
    if 'master_bias' not in master_bias.keywords.keys() and not skip_keywords:
        if 'bias_corrected' in ccd.keywords.keys() and not skip_keywords:
            raise ValueError('This ccd data is already bias subtracted.')
        raise ValueError('master_bias is not a valid master bias frame.')
    else:
        ccd_math('-', ccd, master_bias)
        add_keyword(ccd, 'bias_corrected', 1)
        add_keyword(ccd, 'combine_bias_method', master_bias.keywords['combine_bias_method'])

def combine_flat(**kwargs):
    '''
    Combine flat images in a master_flat, subtracting bias.
    '''
    #TODO: implement with dark correction
    ccd_list = files_or_ccd(**kwargs)
    method = get_method(**kwargs)
    master_bias = kwargs.pop('master_bias')

    if kwargs.get('lacosmic', False):
        lacosmic_subtract(ccd_list)

    for ccd in ccd_list:
        if 'bias_corrected' not in ccd.keywords.keys():
            subtract_bias(ccd, master_bias)

    func = combine_funcs[method]
    result = func(ccd_list)
    add_keyword(result, 'master_flat', 1)
    add_keyword(result, 'combine_flat_method', method)
    ccd_math('/', result, result.data.mean())
    return result

def divide_flat(ccd, master_flat, skip_keywords=False):
    '''
    Divide a CCDData image by a master_flat.
    '''
    if 'master_flat' not in master_flat.keywords.keys() and not skip_keywords:
        if 'flat_corrected' in ccd.keywords.keys() and not skip_keywords:
            raise ValueError('This ccd data is already flat divided.')
        raise ValueError('master_flat is not a valid master flat frame.')
    else:
        ccd_math('/', ccd, master_flat)
        add_keyword(ccd, 'flat_corrected', 1)
        add_keyword(ccd, 'combine_flat_method', master_flat.keywords['combine_flat_method'])
