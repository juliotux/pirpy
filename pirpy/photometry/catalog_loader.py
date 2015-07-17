from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.simbad.core import SimbadObjectIDsResult

from ..log import log

__all__ = ['CatalogLoader', 'order_names']

order_names = ['RMC', 'HD', 'LHA', 'HIP', 'AAVSO', 'ASAS', 'UCAC4', 'TYC', '2MASS', 'OGLE']

class Catalog(object):
    '''
    Stores the configuration of a single catalog.
    '''
    def __init__(self, name,  description='', vizier_table=None, vizier_query=False):
        self._name = name
        self._table = vizier_table

        if vizier_query:
            self._description = Vizier.find_catalogs(name)[vizier_table].description
        else:
            self._description = description

    @property
    def name(self):
        return self._name

    @property
    def description(self):
        return self._description

    def _set_names(self, table, **kwargs):
        raise Exception('Reimplement it, empty catalog.')

    def _keys(self, filter=None):
        raise Exception('Reimplement it, empty catalog.')

    def _load(self, center, radius, filter=None):
        raise Exception('Reimplement it, empty catalog.')

    def get_keys(self, filter=None, **kwargs):
        return self._keys(filter)

    def query(self, center, radius, filter=None, **kwargs):
        return self._load(center, radius, filter)

class QueryData(object):
    '''
    Stores the data of a query.
    '''
    def __init__(self, query, keys):
        self.table = query
        self.keys = keys

################################################################################

class Simbad_Catalog(Catalog):
    def __init__(self):
        Catalog.__init__(self, 'Simbad', description='Simbad online database.')
        self._simbad = None

    def _parse_names_with_priority(self, name, names_list, priority):
        for i in priority:
            for j in names_list:
                if j['ID'][0:len(i)] == i:
                    return j['ID']

        return name

    def _set_names(self, table, **kwargs):
        if self._simbad is None:
            self._simbad = Simbad()

        names_order = kwargs.pop('names_priority', None)
        if names_order is not None:
            names_list = [None]*len(table)

            for i in range(len(names_list)):
                names_list[i] = Simbad.query_objectids_async(table[i]['MAIN_ID'])

            names_list = [self._simbad._parse_result(i, SimbadObjectIDsResult) for i in names_list]

            for i in range(len(table)):
                table[i]['MAIN_ID'] = self._parse_names_with_priority(table[i]['MAIN_ID'], names_list[i], names_order)

        return table

    def _keys(self, filter=None):
        if filter is None:
            log.warn("No filter specified, loading just astrometric data.")
            return {'id_key':'MAIN_ID',
                    'ra_key':'RA',
                    'dec_key':'DEC',}

        return {'id_key':'MAIN_ID',
                'ra_key':'RA',
                'dec_key':'DEC',
                'flux_key':'FLUX_%s' % filter,
                'flux_error_key':'FLUX_ERROR_%s' % filter,
                'flux_bib_key':'FLUX_BIBCODE_%s' % filter,
                'flux_unit_key':'FLUX_UNIT_%s' % filter}

    def _load(self, center, radius, filter=None, names_priority=order_names):
        '''
        Loads the object list from a Simbad query, centered in `center` and
        within a `radius`. It can update existent objects, including names.
        '''
        if self._simbad is None:
            self._simbad = Simbad()

        if filter is None:
            log.warn("No filter specified, loading just astrometric data.")
        else:
            self._simbad.add_votable_fields('fluxdata(%s)' % filter)

        query = self._simbad.query_region(center, radius)
        query = self._set_names(query, names_priority=names_priority)

        if filter is not None:
            self._simbad.remove_votable_fields('fluxdata(%s)' % filter)

        return query

    def query(self, center, radius, filter=None, **kwargs):
        return self._load(center, radius, filter, names_priority=kwargs.pop('names_priority',order_names))

################################################################################
class CatalogLoader(object):
    '''
    This class helps to load a big list of catalogs to photobject colection.
    '''
    def __init__(self, *args, **kwargs):
        self._catalog_list = {'Simbad':Simbad_Catalog()}

    def query_catalog(self, catalog_name, center, radius, filter=None, *args, **kwargs):
        '''
        Load objects from a catalog.
        '''
        if catalog_name in self._catalog_list.keys():
            cat = self._catalog_list[catalog_name]
        else:
            raise ValueError('Catalog named %s not found. Use the list_catalogs() method tu see the available catalogs.' % catalog_name)

        return QueryData(cat.query(center, radius, filter), cat.get_keys())

    def list_catalogs(self):
        '''
        List the available catalogs.
        '''
        for i in self._catalog_list.keys():
            print(i + ' : ' + self._catalog_list[i].description)
    
