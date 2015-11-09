'''
Handle the position and magnitude catalogs.
'''

from astropy.table import Table
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import Angle
from scipy.spatial import cKDTree
import numpy as np

from ..math.list_tools import to_list, match_lengths
from .catalog_loader import CatalogLoader, parse_names_with_priority, resolve_name_from_simbad
from .catalog_loader import order_names as default_order_names
from ..mp import mult_ret, num_threads

from ..log import log

__all__ = ['ObjectCatalog']

coord_type = np.dtype([('id','S32'),('ra','longdouble'),('dec','longdouble'),
                       ('mag','float32'),('mag_err','float32'),('mag_unit','S16')])

check_type = np.dtype([('id','S32'),('ra','longdouble'),('dec','longdouble'),
                       ('mag','float32'),('mag_err','float32'),('mag_unit','S16'),
                       ('cat_id','S32'),('exists',bool)])

def try_float(value):
    '''
    Tries to convert a srting to a float. If not, return false.
    '''
    try:
        float(value)
        return True
    except:
        return False

def coords_from_table(table, ra_key='RA', dec_key='DEC', unit=('hourangle','degree')):
    ra = [None]*len(table)
    dec = [None]*len(table)

    for i in range(len(table)):
        try:
            r = table[i][ra_key]
            d = table[i][dec_key]
        except KeyError:
            raise ValueError("Problems loading tables: wrong keys: ra_key=%s dec_key=%s" %
                             (ra_key, dec_key))

        if not try_float(r) and not try_float(d):
            c = SkyCoord(r, d, frame='icrs', unit=unit)
            ra[i] = c.ra.degree
            dec[i] = c.dec.degree
        else:
            ra[i] = float(r)
            dec[i] = float(d)

    return(ra, dec)

def id_from_table(table, id_key='ID'):
    try:
        id = table[id_key].data
    except KeyError:
        raise ValueError("Problems loading tables: wrong keys: id_key=%s" % id_key)
    return id

def mags_from_table(table, flux_key='FLUX', flux_error_key='FLUX_ERROR', flux_unit_key='FLUX_UNIT'):
    try:
        flux = table[flux_key].data
        erro = table[flux_error_key].data
        unit = table[flux_unit_key].data
    except:
        log.error("Problems loading tables: wrong keys: flux_key=%s flux_error_key=%s flux_unit_key=%s" % (flux_key, flux_error_key, flux_unit_key))
        return [np.nan]*len(table), [np.nan]*len(table), [np.nan]*len(table)
    return flux, erro, unit


class ObjectCatalog(object):
    '''
    Catalog to store and operate with positions and magnitudes to do the photometry.
    '''
    def __init__(self, filter):
        self._cat = np.zeros(0, dtype=coord_type)

        self.filter = filter

        self._kdtree = None
        self._uid_n = 0
        self._catalog_loader = None

    def get_ids(self):
        return self._cat['id']

    def _update_kdtree(self):
        self._kdtree = cKDTree(zip(self._cat['ra'], self._cat['dec']))

    def _new_id(self):
        id = "uid%i" % self._uid_n
        self._uid_n += 1
        return id

    def _id_ra_dec(self, mag_limits=None):
        return self._cat['id'], self._cat['ra'], self._cat['dec']

    def _check_existing(self, id, ra, dec, sep_limit=0.5*u.arcsec):
        exists = (False, 0)
        if id in set(self._cat['id']):
            exists = (True, id)
        elif ra is not None and dec is not None:
            try:
                index, dist = self._kdtree.query((ra, dec), n_jobs=num_threads)
                if dist*u.degree <= Angle(sep_limit):
                    exists = (True, self._cat['id'][index])
            except AttributeError:
                pass
        return exists

    def _update_object(self, obj, update=['name','coordinates','magnitude']):
        if len(update) > 0:
            if 'name' in update and obj['cat_id'] != obj['id']:
                self.rename_star(obj['cat_id'], obj['id'])

            filt = self._cat['id'] == obj['id']
            if 'coordinates' in update:
                if not np.isnan([obj['ra'], obj['dec']]).any():
                    self.set_coordinates(obj['id'], obj['ra'], obj['dec'])

            if 'magnitude' in update:
                if not np.isnan([obj['mag']]).any():
                    self.set_magnitude(obj['id'], obj['mag'], obj['mag_err'], obj['mag_unit'])

    def rename_star(self, old_name, new_name):
        '''
        Renames a star in the list.
        '''
        if new_name is None:
            log.warn('Object %s will not be renamed to None.' % old_name)
        else:
            if not old_name in set(self._cat['id']):
                raise ValueError('Star %s not in the database' % old_name)
            if new_name != old_name:
                for i in range(len(self._cat)):
                    if self._cat['id'][i] == old_name:
                        log.info('Object %s will be renamed to %s.' % (old_name, new_name))
                        self._cat['id'][i] = new_name

    def set_coordinates(self, id, ra, dec):
        '''
        Set the coordinates of a single object.

        Parameters:
            ra : float
                The RA of the object in degrees.
            dec : float
                The DEC of the object in degrees.
        '''
        log.info("Object %s will be updated to the coordinates: ra=%f, dec=%f" % (id, ra, dec))
        for i in range(len(self._cat)):
            if self._cat['id'][i] == id:
                self._cat['ra'][i] = float(ra)
                self._cat['dec'][i] = float(dec)
        
        self._update_kdtree()

    def set_magnitude(self, id, mag, mag_err, mag_unit):
        log.info("Object %s will be updated to catalog magnitude: mag=%f, mag_err=%f, mag_unit=%s" % (str(id), float(mag), float(mag_err), str(mag_unit)))

        for i in range(len(self._cat)):
            if self._cat['id'][i] == id:
                self._cat['mag'][i] = float(mag)
                self._cat['mag_err'][i] = float(mag_err)
                self._cat['mag_unit'][i] = str(mag_unit)

    def get_skycoord(self, id):
        '''
        Returns the SkyCoord instance of the coordinates of an object.
        '''
        for i in self._cat:
            if i['id'] == id:
                return SkyCoord(i['ra'], i['dec'], frame='icrs', unit=('degree','degree'))
        return None

    def get_mag(self, id):
        for i in self._cat:
            if i['id'] == id:
                return i['mag']
        return np.nan

    def add_object(self, id, ra=np.nan, dec=np.nan, mag=np.nan, mag_err=np.nan, mag_unit=None,
                   update=['name','coordinates','magnitude'], sep_limit=0.5*u.arcsec,
                   skip_existence_checking=False):
        '''
        Adds an object to the list.
        '''
        self.add_objects([id], [ra], [dec], [mag], [mag_err], [mag_unit],
                         update=update, sep_limit=sep_limit,
                         skip_existence_checking=skip_existence_checking)

    def add_objects(self, id, ra=[np.nan], dec=[np.nan], mag=[np.nan], mag_err=[np.nan], mag_unit=[None],
                    update=['name','coordinates','magnitude'], sep_limit=0.5*u.arcsec,
                    skip_existence_checking=False):
        id = to_list(id)
        ra = to_list(ra)
        dec = to_list(dec)
        mag = to_list(mag)
        mag_err = to_list(mag_err)
        mag_unit = to_list(mag_unit)

        id, ra, dec, mag, mag_err, mag_unit = match_lengths([id, ra, dec, mag, mag_err, mag_unit], len(id))

        for i in range(len(id)):
            if id[i] == None or id[i] == 'None' or id[i] == '':
                id[i] = self._new_id()

        objs = np.array(zip(id, ra, dec, mag, mag_err, mag_unit, ['']*len(id), [False]*len(id)), dtype=check_type)

        if skip_existence_checking:
            self._cat = np.append(self._cat, np.array(zip(id, ra, dec, mag, mag_err, mag_unit), dtype=coord_type))
        else:
            for i in range(len(objs)):
                tof, id2 = self._check_existing(objs['id'][i], objs['ra'][i], objs['dec'][i], sep_limit=sep_limit)

            self._cat = np.append(self._cat, np.array(zip(objs['id'], objs['ra'], objs['dec'],
                                                          objs['mag'], objs['mag_err'], objs['mag_unit']),
                                                          dtype=coord_type)[~objs['exists']])

            for i in objs[objs['exists']]:
                self._update_object(i, update=update)

        self._update_kdtree()

    def load_objects_from_table(self, table, id_key='ID', ra_key='RA', dec_key='DEC',
                                flux_key='MAG', flux_error_key='MAG_ERR', flux_unit_key='MAG_UNIT', flux_bib_key=None,
                                update=['name','coordinates','magnitude'], format='fixed_width',
                                sep_limit=1*u.arcsec, skip_existence_checking=False):
        '''
        Loads the object list from a table.
        '''
        if not isinstance(table, Table):
            try:
                table = ascii.read(table, format=format)
            except:
                raise ValueError("Problem loading the table. Please give an astropy Table or a filename.")

        id = id_from_table(table, id_key)
        ra, dec = coords_from_table(table, ra_key, dec_key)

        if self.filter is not None:
            mag, err, unit = mags_from_table(table, flux_key, flux_error_key, flux_unit_key)
            self.add_objects(id, ra, dec, mag, err, unit, update=update, sep_limit=sep_limit, skip_existence_checking=skip_existence_checking)
        else:
            self.add_objects(id, ra, dec, update=update, sep_limit=sep_limit, skip_existence_checking=skip_existence_checking)

    def load_objects_from_catalog(self, catalog_name, center, radius,
                                  update=['name','coordinates','magnitude'],
                                  sep_limit=1*u.arcsec, skip_existence_checking=False, **kwargs):
        '''
        Loads the object list from online catalogs.
        '''
        if self._catalog_loader is None:
            self._catalog_loader = CatalogLoader()

        result = self._catalog_loader.query_catalog(catalog_name, center, radius, filter=self.filter, **kwargs)

        self.load_objects_from_table(result.table, update=update, sep_limit=sep_limit, **result.keys)

    def list_online_catalogs(self):
        '''
        List the available catalogs.
        '''
        if self._catalog_loader is None:
            self._catalog_loader = CatalogLoader()

        self._catalog_loader.list_catalogs()

    def match_point(self, ra, dec, r_lim=1*u.arcsec, add_new=False):
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
        ra = to_list(ra)
        dec = to_list(dec)

        if self._kdtree is None:
            for i in zip(ra, dec):
                self.add_object(None, i[0], i[1])

        dist, index = self._kdtree.query(zip(ra, dec), n_jobs=num_threads)

        ids = self._cat['id'][index]
        for i in range(len(ids)):
            if dist[i]*u.degree > Angle(r_lim):
                if add_new:
                    self.add_object(None, ra[i], dec[i])
                    ind, dis = self._kdtree.query((ra[i], dec[i]), n_jobs=num_threads)
                    ids[i] = self._cat[ind]['id']
                else:
                    ids[i] = 'None'
        return ids
