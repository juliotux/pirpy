'''
Implements the photometry storage for a single object.
'''

from __future__ import division

from astropy.table import Table

from ..math.list_tools import to_list, match_lengths

from ..log import log

__all__ = ['PhotObject', 'PhotColection']

class PhotObject(object):
    '''
    The photometry storage for a single object, with fast comparision with
    another.
    '''
    def __init__(self, id, ra, dec, filter = None):
        '''
        Iinitializes a PhotObject py giving an `ra`, `dec` and,
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
    def filter(self):
        return self._filter

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
        return [(j,k) for i,j,k in
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
            log.error('All the arguments must be numbers. Got: jd=%f   sum=%f   error=%f' % (jd, sum, error))
            raise ValueError('Not number values.')

        if filter != self._filter and filter is not None:
            log.error('The filter of the meassurement is different from the filter of the object.')

        if jd in set(self._jds):
            log.warn('There is already a measure for this time.')
        else:
            self._jds.append(jd)
            self._sums.append(sum)
            self._sums_error.append(error)
            log.debug('Added jd=%f   sum=%f   error=%f to object %s' % (jd, sum, error, self._id) )

    def relative_to(self, photobject):
        '''
        Make the relation of the photometry of this object to other.
        '''
        fjd = []
        fphot = []
        ferror = []

        for jd1, pht1, err1 in zip(self._jds, self._sums, self._sums_error):
            if jd1 in set(photobject.jd):
                pht2, err2 = photobject.get_sum(jd1)
                pht = pht1/pht2
                err = pht*((err1/pht1) + (err2/pht2))
                fjd.append(jd1)
                fphot.append(pht)
                ferror.append(err)

        indx = np.argsort(fjd)
        return fjd[indx], fphot[indx], ferror[indx]


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
        id = []
        ra = []
        dec = []

        for i in self._list.keys():
            id.append(i)
            ra.append(self._list[i].ra)
            dec.append(self._list[i].dec)
        return Table([id, ra, dec], names=('ID', 'RA', 'DEC'))

    @property
    def get_photometry(self):
        id = []
        jd = []
        flux = []
        error = []

        for i in self._list.keys():
            for j, f, e in zip(self._list[i].jd,
                               self._list[i].sums,
                               self._list[i].sums_error):
                id.append(i)
                jd.append(j)
                flux.append(f)
                error.append(e)

        return Table([id, jd, flux, error], names=('ID','JD','FLUX','FLUX_ERROR'))

    def get_object(self, id):
        '''
        Returns a complete object with an `id`.
        '''
        return self._list[id]

    def add_objects(self, id, ra=None, dec=None):
        '''
        Adds an object to the list.
        '''
        if id in self._list.keys():
            raise ValueError('The object %s is already in the colection' % id)
        else:
            self._list[id] = PhotObject(id, ra, dec, self._filter)
            log.debug('Added the object %s to PhotColection' % id)

    def add_results(self, jd, id, flux, flux_error, ra=None, dec=None, add_if_not_exists=True):
        '''
        Append a result list to the objects, guided by id.
        '''
        jd = to_list(jd)
        id = to_list(id)
        ra = to_list(ra)
        dec = to_list(dec)
        flux = to_list(flux)
        flux_error = to_list(flux_error)

        jd, id, ra, dec, flux, flux_error = match_lengths([jd, id, ra, dec, flux, flux_error], len(id))

        for i in range(len(id)):
            if id[i] not in self.ids and add_if_not_exists:
                self.add_objects(id[i], ra[i], dec[i])
            if id[i] in self.ids:
                self._list[id[i]].add_result(jd[i], flux[i], flux_error[i])

    def compare_objects(self, id1, id2):
        '''
        Compare the photometry of 2 objects.
        '''
        return self._list[id1].relative_to(self._list[id2])

    def load_from_table(self, table, id_key='ID', jd_key='JD', flux_key='FLUX', error_key='FLUX_ERROR'):
        '''
        Loads the photometry from a ~astropy.table~
        '''
        if not isinstance(table, Table):
            from astropy.io import ascii
            try:
                table = ascii.read(table)
            except:
                raise ValueError("Problem loading the table. Please give an astropy.table.Table instance or a filename.")

        for i in table:
            try:
                id = i[id_key]
                jd = i[jd_key]
                fl = i[flux_key]
                er = i[error_key]
            except KeyError:
                raise ValueError("Problems loading tables: wrong keys: id_key=%s jd_key=%s flux_key=%s error_key=%s" %
                                 (id_key, jd_key, flux_key, error_key))

    def write_to(self, fname):
        '''
        Save the results to a file.
        '''
        from astropy.io import ascii
        result = self.get_photometry
        ascii.write(result, fname)
