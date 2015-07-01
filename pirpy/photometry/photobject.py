'''
Implements the photometry storage for a single object.
'''

from __future__ import division

from astropy.table import Table

from ..math.list_tools import to_list, match_lengths

__all__ = ['PhotObject', 'PhotColection']

class PhotObject(object):
    '''
    The photometry storage for a single object, with fast comparision with
    another.
    '''
    def __init__(self, ra, dec, filter = None):
        '''
        Iinitializes a PhotObject py giving an `ra`, `dec` and,
        optionalliy, a `filter`.
        '''
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
        return self._jds

    def get_sum(self, jd):
        return [j,k for i,j,k in
                zip(self._jds, self._sums, self._sums_error)
                if i == jd][0]

    def add_result(self, jd, sum, error=0.0, filter=None):
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
        Make the relation of the photometry of this object to other.
        '''
        #TODO: Make the algorithm faster
        jd = phot = error = []

        for i in range(len(self._jds)):
            jd1 = self._jds[1]
            if jd1 in set(photobject.jd):
                phot1 = self._sums[i]
                error1 = self._sums_error[i]
                phot2, error2 = photobject.get_sum(jd1)
                photf = phot1/phot2
                errorf = photf*((error1/phot1) + (error2/phot2))
                jd.append(jd1)
                phot.append(photf)
                error.append(errorf)

        return jd, phot, error


class PhotColection(object):
    '''
    A colection to store a lot of photobjects.
    '''
    def __init__(self, filter=None):
        self._list = {}
        self._filter = None

    @property
    def ids(self):
        return self._list.keys()

    @property
    def objects(self):
        id = ra = dec = []
        for i in self._list.keys():
            id.append(i)
            ra.append(self._list[i].ra)
            dec.append(self._list[i].dec)
        return Table([id, ra, dec],
                     names=['ID', 'RA', 'DEC'],
                     unit=[None, 'degree', 'degree'])

    @property
    def get_photometry(self):
        id = jd = flux = flux_error = []

        for i in self._list.keys():
            for j, f, e in zip(self._list[i].jd,
                               self._list[i].sums,
                               self._list[i].sums_error):
                id.append(i)
                jd.append(j)
                flux.append(f)
                flux_error.append(e)

        return Table([id, jd, flux, flux_error],
                     names = ('ID', 'JD', 'FLUX', 'FLUX_ERROR'))

    def add_objects(self, id, ra=None, dec=None):
        '''
        Adds a object to the list.
        '''
        if id in self._list.keys():
            raise ValueError('The object ' + str(id) + ' is already in the colection')
        else:
            self._list[id] = PhotObject(ra, dec, self._filter)

    def add_results(self, jd, id, ra, dec, flux, flux_error, add_if_not_exists=True):
        '''
        Append a result list to the objects, guided by id.
        '''
        jd = to_list(jd)
        ra = to_list(ra)
        dec = to_list(dec)
        flux = to_list(flux)
        flux_error = to_list(flux_error)
        flux_flag = to_list(flux_flag)

        jd, id, ra, dec, flux, flux_error = match_lengths([jd, id, ra, dec, flux, flux_error], len(id))

        for i in range(len(id)):
            if id[i] not in self.ids and add_if_not_exists:
                self.add_objects(id[i], ra[i], dec[i])
            if id[i] in self.ids:
                self._list[i].add_result(jd[i], flux[i], flux_error[i])

    def compare_objects(self, id1, id2):
        '''
        Compare the photometry of 2 objects.
        '''
        return self._list[id1].relative_to(self._list[id2])
