'''
Implements the photometry storage for a single object.
'''

from __future__ import division

from astropy.table import Table
from astropy.io import ascii
from astropy import units as u
from scipy.spatial import cKDTree

from ..math.list_tools import to_list, match_lengths

from ..log import log

__all__ = ['PhotObject', 'PhotColection']

def get_simbad_config(filter):
    '''
    Get the configuration names for the Simbad tables.
    '''
    return {'id_key':'MAIN_ID',
            'ra_key':'RA',
            'dec_key':'DEC',
            'flux_key':'FLUX_'+filter,
            'flux_error_key':'FLUX_ERROR_'+filter,
            'flux_bib_key':'FLUX_BIBCODE_'+filter,
            'flux_unit_key':'FLUX_UNIT_'+filter}

class PhotObject(object):
    '''
    The photometry storage for a single object, with fast comparision with
    another.
    '''
    def __init__(self, id, ra, dec, filter = None,
                 cat_mag=None, cat_mag_err=None, cat_mag_unit=None):
        '''
        Iinitializes a PhotObject py giving an `ra`, `dec` and,
        optionalliy, a `filter`.
        '''
        self._id = id #Id
        self._ra = ra #RA
        self._dec = dec #Dec
        self._filter = filter #Filter

        self._sums = [] #Meassured flux
        self._sums_error = [] #Error of each meassurement
        self._jds = [] #Julian dates of the meassures

        self._cat_mag = cat_mag #Catalog magnitude
        self._cat_err = cat_mag_err #Error of the magnitude
        self._cat_unit = cat_mag_unit #Unit of the magnitude

    @property
    def ra(self):
        return self._ra

    @property
    def dec(self):
        return self._dec

    @property
    def id(self):
        return self.id

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

    @property
    def cat_mag(self):
        return self._cat_mag

    @proprty
    def cat_err(self):
        return self._cat_err

    @property
    def cat_unit(self):
        return self._cat_unit

    def set_ra(self, ra):
        self._ra = ra

    def set_dec(self, dec):
        self._dec = dec

    def set_id(self, id):
        self._id = id

    def set_mag(self, mag, mag_err=None, mag_unit=None):
        self._cat_mag = mag
        self._cat_err = mag_err
        self._cat_unit = mag_unit

    def set_catalog_magnitudes(self, mag, err=None, bib=None):
        self._cat_mag = mag
        if err is not None:
            self._cat_err = err
        if bib is not None:
            self._cat_bib = bib

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
            log.warn('There is already a measure for this time. id=%s, jd=%f' % (self._id, jd))
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
        self._filter = filter

        self._kdtree = None
        self._kd_ids = None
        self._mag_bib_list = set([])

        self._uid_n = 0

    @property
    def ids(self):
        return self._list.keys()

    def _id_ra_dec(self):
        id = []
        ra = []
        dec = []

        for i in self._list.keys():
            id.append(i)
            ra.append(self._list[i].ra)
            dec.append(self._list[i].dec)

        return id, ra, dec

    def _update_kdtree(self):
        self._kd_ids, ra, dec = _id_ra_dec()
        self._kdtree = cKDTree(zip(ra, dec))

    def _query_nn(self, ra, dec):
        '''
        Query the nearest neighbors of a point with coordinates ra and dec.
        '''
        dist, index = self._kdtree.query([ra, dec])
        log.debug("The nearest neighbor of the %f,%f point is the %i: %f, %f, with a distance of %f"
                  % (ra, dec, index, self._ra[index], self._dec[index], dist))
        return(index, dist)

    def _new_id(self):
        id = "uid%i" % self._uid_n
        self._uid_n += 1
        return id

    @property
    def objects(self):
        id, ra, dec = self._id_ra_dec()
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

    def set_filter(self, filter):
        '''
        Set the curent filter.
        '''
        self._filter = filter

    def rename_star(self, old_name, new_name):
        '''
        Renames a star in the list.
        '''
        self._list[new_name] = self._list.pop(old_name)
        self._list[new_name].set_id(new_name)

    def add_object(self, id, ra=None, dec=None, mag=None, mag_err=None, mag_unit=None,
                   update_if_exists = True, update_names=True, sep_limit='00d00m01s'):
        '''
        Adds an object to the list.
        '''

        exists = False

        if id is None:
            id = self._new_id()

        id2 = id

        if id in self._list.keys():
            exists = True
        elif ra is not None and dec is not None:
            index, dist = self._query_nn([ra, dec])
            if dist <= 1.0/3600.0:
                exists = True
                id2 = self._kd_ids[index]

        if not exists:
            self._list[id] = PhotObject(id, ra, dec, self._filter, mag, mag_err, mag_unit)
            log.debug('Added the object %s to PhotColection' % id)
        else:
            if update_if_exists:
                if update_names and id2 != id:
                    self.rename_star(id2, id)
                else:
                    id = id2

                if (ra is not None and dec is not None) and
                   ((self._list[id].ra is None or self._list[id].dec is None) or
                    update_if_exists):
                    log.info("Object %s will be updated to the coordinates: ra=%f, dec=%f" % (id, ra, dec))
                    self._list[id].set_ra(ra)
                    self._list[id].set_dec(dec)
                if (mag is not None and mag_err is not None and mag_unit is not None) and update_if_exists:
                    log.info("Object %s will be updated to catalog magnitude: mag=%f, mag_err=%f, mag_unit=%f")
                    self._list[id].set_mag(mag, mag_err, mag_unit)
            else:
                log.error("Object %s already exists and not updated" % id)

        self._update_kdtree()

    def get_object(self, id):
        '''
        Returns a complete object with an `id`.
        '''
        return self._list[id]

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
                self.add_object(id[i], ra[i], dec[i])
            if id[i] in self.ids:
                self._list[id[i]].add_result(jd[i], flux[i], flux_error[i])

    def relative_light_curve(self, id1, id2):
        '''
        Compare the photometry of 2 objects.
        '''
        return self._list[id1].relative_to(self._list[id2])

    def load_objects_from_table(self, table, id_key='ID', ra_key='RA', dec_key='DEC',
                                flux_key=None, flux_error_key=None, flux_unit_key=None, flux_bib_key=None,
                                update_if_exists = True, update_names=True,
                                sep_limit='00d00m01s'):
        '''
        Loads the object list from a table.
        '''
        if not isinstance(table, Table):
            try:
                table = ascii.read(table)
            except:
                raise ValueError("Problem loading the table. Please give an astropy Table or a filename.")

        for i in table:
            try:
                id = i[id_key]
                ra = i[ra_key]
                dec = i[dec_key]
            except KeyError:
                raise ValueError("Problems loading tables: wrong keys: id_key=%s ra_key=%s dec_key=%s" %
                                 (id_key, ra_key, dec_key))

            if (flux_key is not None) and (flux_error_key is not None) and (flux_unit_key is not None):
                try:
                    mag = i[flux_key]
                    err = i[flux_error_key]
                    unit = i[flux_unit_key]
                except KeyError:
                    raise ValueError("Problems loading tables: wrong keys: flux_key=%s flux_error_key=%s flux_bib_key=%s flux_unit_key=%s" %
                                     (flux_key, flux_error_key, flux_unit_key))
                self.add_object(id, ra, dec, mag, err, unit)
                if flux_bib_key is not None:
                    try:
                        self._mag_bib_list.add(i[flux_bib_key])
                    except:
                        pass
            else:
                self.add_object(id, ra, dec)

    def load_objects_from_simbad(self, center, radius, update_if_exists=True, update_names=True, sep_limit='00d00m01s'):
        '''
        Loads the object list from a Simbad query, centered in `center` and
        within a `radius`. It can update existent objects, including names.

        The source identification is did using the `sep_limit`.
        '''
        from astroquery.simbad import Simbad

        s = Simbad()

        if self._filter is not None:
            s.add_votable_field('fluxdata(%s)' % self._filter)

        query = s.query_region(center, radius)

        self.load_objects_from_table(query, **get_simbad_config(self._filter),
                                     update_if_exists=update_if_exists,
                                     update_names=update_names,
                                     sep_limit=sep_limit)

    def load_photometry_from_table(self, table, id_key='ID', jd_key='JD', flux_key='FLUX', error_key='FLUX_ERROR'):
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
            self.add_results(jd, id, fl, er, add_if_not_exists=True)

    def save_objects_catalog(self, fname):
        '''
        Save the objects catalog to a file.
        '''
        ascii.write(self.objects, fname_objects)

    def save_photometry(self, fname):
        '''
        Save the results to a file.
        '''
        ascii.write(self.get_photometry, fname)    def query_nn(self, ra, dec):
        '''
        Query the nearest neighbors of a point with coordinates ra and dec.
        '''
        dist, index = self._kdtree.query([ra, dec])
        log.debug("The nearest neighbor of the %f,%f point is the %i: %f, %f, with a distance of %f"
                  % (ra, dec, index, self._ra[index], self._dec[index], dist))
        return(index, dist)

    def match_point(self, ra, dec, r_lim='00d00m01s', add_new=False):
        '''
        Match a point with coordinates ra and dec in your catalog. This point
        must be inside the r_lim angle. Otherwise, returns None, or add a new uid
        point if specified.

        Parameters:
            ra : float
                The RA coordinate of the point.
            dec : float
                The DEC coordinate of the point.
            r_lim : ~astropy.quantity~
                The limmit angle to the match. Default: 2 arcsec.
            add_new : bool
                If true, the code will add a new point in the catalog with the
                unmatched coordinates.

        Returns:
            string :
                The ID of the matched object. Returns None if no object matched
                the point.
        '''
        if self._kdtree is None:
            self.add_objects(None, ra, dec)

        index, dist = self._query_nn(ra, dec)

        if dist*u.degree <= r_lim:
            return self._kd_ids[index]
        else:
            if add_new:
                self.add_objects(None, ra, dec)
                return self._kd_ids[-1]
            else:
                return None
