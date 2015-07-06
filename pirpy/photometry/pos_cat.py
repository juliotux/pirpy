from scipy.spatial import cKDTree
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import astropy.units as u
import numpy as np

from ..math.list_tools import to_list

from ..log import log

'''
Module to match a point with RA,DEC coordinates into a catalog. If not match
inside a distance R, it can create a new object in the catalog.
'''

__all__ = ['PositionCatalog']

class PositionCatalog(object):
    '''
    High level object to handle position catalogs to use in photometry.
    '''
    def __init__(self,**kwargs):
        self._ra = []
        self._dec = []
        self._id = []

        self._uid_n = 0

        self._args = list(kwargs.keys())
        if 'table' in self._args:
            self._coords_and_id_from_simbad_table(kwargs['table'])

        self._kdtree = None

    @property
    def ra(self):
        return self._ra

    @property
    def dec(self):
        return self._dec

    @property
    def id(self):
        return self._id

    def add_objects(self, id, ra, dec, pass_if_exists=True):
        '''
        Add objects to the PositionCatalog.

        Parameters:
            id : string, list or array
                The identifications of the objects.
            ra : float or array
                The RA coordinates of the objects.
            dec : float or array
                The DEC coordinates of the objects.
        '''
        id = to_list(id)
        ra = to_list(ra)
        dec = to_list(dec)

        if len(ra) != len(dec) or len(ra) != len(id):
            raise ValueError('The id, ra and dec variables must have the same dimentions')

        #check if the user is trying to add an existent object.
        for i in range(len(id)):
            if id[i] not in set(self._id):
                id2 = ''
                if id[i] is not None:
                    id2 = id[i]
                else:
                    id2 = 'uid' + str(self._uid_n)
                    self._uid_n += 1
                self._id.append(id2)
                self._ra.append(float(ra[i]))
                self._dec.append(float(dec[i]))
                log.debug("Added the object %s: %f,%f to the position catalog"
                          % (id2, float(ra[i]), float(dec[i])))
            else:
                if not pass_if_exists:
                    raise ValueError("The object %s is already in this position catalog." % (id[i]))

        self._kdtree = cKDTree(np.dstack([self._ra, self._dec])[0])

    def _try_float(self, value):
        '''
        Tries to convert a srting to a float. If not, return false.
        '''
        try:
            float(value)
            return True
        except:
            return False

    def _coords_from_table(self, table, ra_key='RA', dec_key='DEC'):
        for i in range(len(table)):
            try:
                ra = table[i][ra_key]
                dec = table[i][dec_key]
            except:
                raise ValueError('The column names for RA and DEC are wrong. ' +
                                 'The values gived by user are: '
                                 + str(ra_key) + ' ' + str(dec_key))
            if not self._try_float(ra) and not self._try_float(dec):
                c = SkyCoord(ra, dec, frame='icrs', unit=('hourangle','degree'))
                table[i][ra_key] = c.ra.degree
                table[i][dec_key] = c.dec.degree

        return(table[ra_key].data.data, table[dec_key].data.data)

    def _coords_and_id_from_table(self, table, id_key='MAIN_ID', ra_key='RA', dec_key='DEC'):
        '''
        Returns arrays with id, ra and dec from an astropy.table in a format like
        the queried from Simbad via ~astroquery~.
        '''
        try:
            id = table[id_key].data.data
        except:
            raise KeyError('The column name for id is wrong. The value gived by user is: %s' % id_key)

        ra, dec = self._coords_from_simbad_table(table, ra_key, dec_key)
        return id, ra, dec

    def add_from_table(self, table, id_key='MAIN_ID', ra_key='RA', dec_key='DEC',
                       pass_if_exists=True):
        '''
        Loads a ~astropy.table~ and append it to the actual catalog.
        '''
        id, ra, dec = self._coords_and_id_from_table(table, id_key, ra_key, dec_key)
        self.add_objects(id, ra, dec, pass_if_exists=pass_if_exists)

    def query_nn(self, ra, dec):
        '''
        Query the nearest neighbors of a point with coordinates ra and dec.
        '''
        dist, index = self._kdtree.query([ra, dec])
        log.debug("The nearest neighbor of the %f,%f point is the %i: %f, %f, with a distance of %f"
                  % (ra, dec, index, self._ra[index], self._dec[index], dist))
        return(index, dist)

    def match_point(self, ra, dec, r_lim=2*u.arcsec, add_new=False):
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

        index, dist = self.query_nn(ra, dec)

        if dist*u.degree <= r_lim:
            log.debug("Matched %f,%f point with %s: %f, %f"
                      % (ra, dec, self._id[index], self._ra[index], self._dec[index]))
            return self._id[index]
        else:
            if add_new:
                log.debug("Not matched, addind new object.")
                self.add_objects(None, ra, dec)
                return self._id[-1]
            else:
                return None

    def add_from_simbad(self, center, radius, pass_if_exists=True):
        '''
        Loads a catalog from an ~astroquery.simbad~ query.

        Parameters:
            center : string or ~astropy.coordinates.SkyCoord~
                The center of the field to query.
            radius : string or ~astropy.units~
                The maximum angle to query objects.
        '''
        from astroquery.simbad import Simbad

        self.add_from_table(Simbad.query_region(center, radius),
                            pass_if_exists=pass_if_exists)

    def write_to(self, fname):
        '''
        Write the actual catalog to a file.
        '''
        table = Table([self._id, self._ra, self._dec],
                      names=('MAIN_ID', 'RA', 'DEC'))

        ascii.write(table, fname)

