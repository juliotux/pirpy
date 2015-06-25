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
        Iinitializes a PhotObject py giving an `id`, `ra`, `dec` and,
        optionalliy, a `filter`.
        '''
        self._id = id
        self._ra = ra
        self._dec = dec
        self._filter = filter

        self._sums = []
        self._sums_error = []
        self._jds = []

    @property
    def ra(self):
        return self._ra

    @property
    def dec(self):
        return self._dec

    @property
    def id(self):
        return self._id

    @property
    def filter(self):
        return self.filter

    @property
    def sums(self):
        return self._sums

    @property
    def sums_error(self):
        return self._sums_error

    @property
    def jd(self):
        return self._jd

    def append_result(self, jd, sum, error=None, filter=None):
        '''
        Appends a result [jd, sum, error] to the object.

        Parameters:
            jd : float
                The julian date or another type of time marking.
            sum : float
                The value of meassured photometry.
            error : float
                The error in the photometry meassurement.
        '''
        try:
            jd = float(jd)
            sum = float(sum)
            error = float(error)
        except:
            raise ValueError('All the arguments must be numbers.')

        if filter != self._filter and filter is not None:
            raise ValueError('The filter of the meassurement is different from' +
                             ' the filter of the object').

        if jd in set(self._jds):
            raise ValueError('There is already a measure for this time.')

        self._jds.append(jd)
        self._sums.append(sum)
        self._sums_error.append(error)

    def relative_to(self, photobject):
        '''
        TODO
        '''
