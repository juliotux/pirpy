'''
Arithmetic operations with CCDData.
'''

from ccdproc import CCDData
from ..mp import mult_ret

def __mutiply(ccddata):
    # [0] : CCDData
    # [1] : value to multiply
    return ccddata[0].multiply(ccddata[1])

def __divide(ccddata):
    # [0] : CCDData
    # [1] : value to divide
    return ccddata[0].divide(ccddata[1])

def __subtract(ccddata):
    # [0] : CCDData
    # [1] : value to subtract
    return ccddata[0].subtract(ccddata[1])

def __add_ccd(ccddata):
    # [0] : CCDData
    # [1] : value to sum
    return ccddata[0].add(ccddata[1])

def op(op, ccdlist, value, **kwargs):
    '''
    Executes the mathematical `op` operation with `ccdlist` and `value`.

    Parameters:
        op : string
            It may assume the values:
                - ['subtract','minus','-'] : Subtract operation.
                - ['add','plus','+'] : Sum operation.
                - ['multiply','*'] : Multiply operation.
                - ['divide','/'] : Divide operation.
        ccdlist : list of ccdproc.CCDData
            The list containing the data to multiply.
        value : float or ccdproc.CCDData
            The value, or array with values, to multiply the data.
        nprocess (optional) : int
            Maximum number of parallel processes. Default: 1

    Returns:
        list of ccdproc.CCDData:
            A list with the multiplied data.

    '''
    ccdlist = verify_ccdlist(ccdlist)
    for i in xrange(len(ccdlist)):
        ccdlist[i] = [ccdlist[i], value]

    if op is '-' or op is 'minus' or op is 'suptract':
        return mult_ret(__subtract, ccdlist, kwargs.pop('nprocess',1))
    if op is '+' or op is 'plus' or op is 'add':
        return mult_ret(__add, ccdlist, kwargs.pop('nprocess',1))
    if op is '*' or op is 'multiply':
        return mult_ret(__multiply, ccdlist, kwargs.pop('nprocess',1))
    if op is '/' or op is 'divide':
        return mult_ret(__divide, ccdlist, kwargs.pop('nprocess',1))

    raise ValueError('Unrecognized ' + str(op) + ' operation.')
