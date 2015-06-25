'''
Implements the photometry storage for a single object.
'''

class PhotObject(object):
    '''
    The photometry storage for a single object, with fast comparision with
    another.
    '''
    def __init__(self, id, ra, dec, filter = None):
        '''
        TODO
        '''
        self._id = id
        self._ra = ra
        self._dec = dec
        self._filter = filter

        self._sums = []
        self._sums_error - []
        self._jds = []

    def append_result(self, jd, sum, error=None):
        '''
        TODO
        '''

    def relative_to(self, photobject):
        '''
        TODO
        '''
