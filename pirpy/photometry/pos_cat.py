from scipy.spatial import cKDTree
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import astropy.units as u
import numpy as np

'''
Module to match a point with RA,DEC coordinates into a catalog. If not match
inside a distance R, it can create a new object in the catalog.
'''

class PositionCatalog(object):
    '''
    High level object to handle position catalogs to use in photometry.
    '''
    #TODO: modificar as funcoes de carregamento de tabelas para possibilitar
    #      persinalização das chaves usadas.
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

    def add_objects(self, id, ra, dec):
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
        if not isinstance(id, (tuple, list, set, np.ndarray)):
            id = [id]
        if not isinstance(ra, (tuple, list, set, np.ndarray)):
            ra = [ra]
        if not isinstance(dec, (tuple, list, set, np.ndarray)):
            dec = [dec]

        if len(ra) != len(dec) or len(ra) != len(id):
            raise ValueError('The id, ra and dec variables must have the same dimentions')

        #check if the user is trying to add an existent object. If yes, skip it
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

        self._kdtree = cKDTree(np.dstack([ra, dec])[0])

    def _try_float(self, value):
        '''
        Tries to convert a srting to a float. If not, return false.
        '''
        try:
            float(value)
            return True
        except:
            return False

    def _coords_from_simbad_table(self, table):
        for i in range(len(table)):
            ra = table[i]['RA']
            dec = table[i]['DEC']
            if not self._try_float(ra) and not self._try_float(dec):
                c = SkyCoord(ra, dec, frame='icrs', unit=('hourangle','degree'))
                table[i]['RA'] = c.ra.degree
                table[i]['DEC'] = c.dec.degree
        return(table['RA'].data.data, table['DEC'].data.data)

    def _coords_and_id_from_simbad_table(self, table):
        '''
        Returns arrays with id, ra and dec from an astropy.table in a format like
        the queried from Simbad via ~astroquery~.
        '''
        try:
            id = table['MAIN_ID'].data.data
            ra, dec = self._coords_from_simbad_table(table)
        except:
            raise KeyError('The astropy.table column names to read are' +
                           ' MAIN_ID, RA and DEC. Please check your table.')

        self.add_objects(id, ra, dec)

    def query_nn(self, ra, dec):
        '''
        Query the nearest neighbors of a point with coordinates ra and dec.
        '''
        dist, index = self._kdtree.query([ra, dec])
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
        index, dist = self.query_nn(ra, dec)

        if dist*u.degree <= r_lim:
            return self._id[index]
        else:
            if add_new:
                add_objects(None, ra, dec)
                return self._id[-1]
            else:
                return None

    def from_simbad(self, center, radius):
        '''
        Loads a catalog from an ~astroquery.simbad~ query.

        Parameters:
            center : string or ~astropy.coordinates.SkyCoord~
                The center of the field to query.
            radius : string or ~astropy.units~
                The maximum angle to query objects.
        '''
        from astroquery.simbad import Simbad

        query = Simbad.query_region(center, radius)
        self._coords_and_id_from_simbad_table(query)

    def write_to(self, fname):
        '''
        Write the actual catalog to a file.
        '''
        table = Table([self._id, self._ra, self._dec],
                      names=('MAIN_ID', 'RA', 'DEC'))

        ascii.write(table, fname)

class PhotometryCatalog(object):
    '''
    High level object to handle catalogs of photometric standart stars to use in
    the comparision in relative photometry.
    '''
    def __init__(self, filter, *args, **kwargs):
        self._id = []
        self._flux = []
        self._flux_error = []
        self._flux_unit = []
        self._flux_bib = []
        #The filter argument forces the catalog to be in a single filter only
        self._filter = filter

    def _mag2flux(self, mag, offset=1.0):
        '''
        Converts a magnitude to a linear flux.
        '''
        return offset * 10**(-0.4*mag)

    def _flux2mag(self, flux, offset = 0.0):
        '''
        Converts a linear flux to magnitude.
        '''
        return -2.5*np.log10(flux) + offset

    def add_objects(self, id, flux, unit, flux_error = None, bibcode = None):
        '''
        Add objects to the catalog.

        Parameters:
            id : string, list or array
                The identifications of the objects.
            flux : float or array
                The fluxes of the objects.
            flux_error : float or array
                The error of the flux.
            unit : float or array
                The unit of the fluxes.
        '''
        if not isinstance(id, (tuple, list, set, np.ndarray)):
            id = [id]
        if not isinstance(flux, (tuple, list, set, np.ndarray)):
            flux = [flux]
        if not isinstance(flux_error, (tuple, list, set, np.ndarray)):
            flux_error = [flux_error]
        if not isinstance(unit, (tuple, list, set, np.ndarray)):
            unit = [unit]*len(flux)
        if not isinstance(bibcode, (tuple, list, set, np.ndarray)):
            bibcode = [bibcode]*len(flux)

        if len(id) != len(flux):
            raise ValueError('The id flux variables must have the same dimentions')

        #check if the user is trying to add an existent object. If yes, skip it
        for i in range(len(id)):
            if id[i] not in set(self._id):
                id2 = ''
                if id[i] is not None:
                    id2 = id[i]
                else:
                    id2 = 'uid' + str(self._uid_n)
                    self._uid_n += 1
                self._id.append(id2)
                self._flux.append(float(flux[i]))
                self._flux_error.append(float(flux_error[i]))
                self._unit.append(unit[i])
                self._flux_bib.append(bibcode[i])

    def _flux_and_id_from_table(self, table, id_key, flux_key, unit_key,
                                error_key=None, bib_key=None):
        '''
        Add objects from a ~astropy.table~.

        Parameters:
            table : ~astropy.table.Table~
                The table with the data.
            id_key : string
                The name of the column that contains the ID of the objects.
            flux_key : sting
                The name of the column that contains the flux of the objects.
            unit_key : string
                The name of the column that contains the flux unit of the object.
            error_key : string (optional)
                The name of the column that contains the flux error of the object.
            bib_key : string (optional)
                The name of the column that contains the bibcode of the measure.
        '''
        id = []
        flux = []
        flux_err = []
        unit = []
        bib = []
        #TODO: continuar aqui


    def from_simbad(self, center, radius):
        '''
        Add objects from an ~astroquery.simbad~ query.
        '''
        from astroquery.simbad import Simbad

        s = Simbad()
        s.add_votable_fields('fluxdata(' + str(self._filter) + ')')

        self._flux_and_id_from_simbad_table(s.query_region(center, radius))

    def generate_bibcode_list(self, fname=None):
        '''
        Generates a bibcode list of the references and possible write it to a file.
        '''
        biblist = set(self._flux_bib) - set([None, '', ' ', '-', '--', '---'])
        #TODO: Continuar aqui

    def compare_object(self, counts, refs_id, refs_counts,
                       counts_error=None, ref_counts_error=None):
        '''
        Compare a object with a reference to get the diferencial photometry.
        '''
        
