from ccdproc import CCDData

def verify_ccdlist(ccdlist):
    '''
    Verify the ccdlist variables if it is a single ccdproc.CCDData instance or
    if is a list of this type of data. Returns a list with the CCDData
    '''
    if isinstance(ccdlist, list):
        for i in ccdlist:
            if not isinstance(i, CCDData):
                raise ValueError('You must give a list of ccdproc.CCDData' +
                                 ' instances in ccdlist variable.')
        return ccdlist

    if isinstance(ccdlist, CCDData):
        return list([ccdlist])

