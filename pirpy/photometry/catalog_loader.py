from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.irsa import Irsa

from astropy.coordinates.angles import Angle
from astropy.coordinates import SkyCoord
from astropy.table import Column

from ..mp import mult_ret

from ..log import log

__all__ = ['CatalogLoader', 'order_names']

order_names = ['RMC', 'HD', 'LHA', 'HIP', 'AAVSO', 'ASAS', 'UCAC4', 'TYC']

Vizier.columns = ['all']

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
        return self._keys(filter, **kwargs)

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

def simbad_query_ids(name):
    return Simbad.query_objectids(name)

def simbad_query_object0(center):
    return Simbad.query_region(center)[0]

def mult_simbad_query_ids(names, nprocess=10):
    return mult_ret(simbad_query_ids, names, nprocess=nprocess)

def mult_simbad_query_object0(centers, nprocess=10):
    return mult_ret(simbad_query_object0, centers, nprocess=nprocess)

def parse_names_with_priority(name, names_list=None, names_priority=order_names):
    if names_list is None:
        return name

    for i in names_priority:
        for j in names_list:
            if j['ID'][0:len(i)] == i:
                return j['ID']

    return name

def resolve_name_from_simbad(sipa):
    '''
    Resolves a name of one star from Simbad, following a priority order.
    '''
    skycoord = sipa[0]
    id = sipa[1]
    priority = sipa[2]
    angle_lim = sipa[3]

    q = Simbad.query_region(skycoord)[0]
    if skycoord.separation(SkyCoord(q['RA'], q['DEC'], frame='icrs', unit=('hourangle','degree'))) <= Angle(angle_lim):
        names_list = Simbad.query_objectids(q['MAIN_ID'])
        return id, parse_names_with_priority(id, names_list, names_priority=priority)
    else:
        return id, None

class Simbad_Catalog(Catalog):
    def __init__(self):
        Catalog.__init__(self, 'Simbad', description='Simbad online database.')
        self._simbad = None

    def _set_names(self, table, **kwargs):
        names_order = kwargs.pop('names_priority', None)
        if names_order is not None:
            names_list = mult_simbad_query_ids(table['MAIN_ID'].data.tolist())

            for i in range(len(table)):
                table[i]['MAIN_ID'] = parse_names_with_priority(table[i]['MAIN_ID'], names_list[i], names_order)

        return table

    def _keys(self, filter=None):
        if filter is None:
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

    def _load(self, center, radius, filter=None, names_priority=None):
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

class UCAC4_Catalog(Catalog):
    def __init__(self):
        Catalog.__init__(self, 'UCAC4', vizier_table='I/322A', vizier_query=True)
        Vizier.ROW_LIMIT = -1

    def _set_names(self, table, **kwargs):
        return  ["UCAC4 %s" % i['UCAC4'] for i in table]

    def _keys(self, filter=None):
        if filter is None:
            log.warn("No filter specified, loading just astrometric data.")
            return {'id_key':'ID',
                    'ra_key':'RAJ2000',
                    'dec_key':'DEJ2000'}

        if filter not in {'V','B','H','J','K','g','r','i'}:
            raise ValueError("UCAC4 catalog only have the following filters: B V J H K g r i")

        return {'id_key':'ID',
                'ra_key':'RAJ2000',
                'dec_key':'DEJ2000',
                'flux_key':'%smag' % filter,
                'flux_error_key':'e_%smag' % filter,
                'flux_bib_key':'bib',
                'flux_unit_key':'unit'}

    def _load(self, center, radius, filter=None):
        table = Vizier.query_region(center, radius=Angle(radius), catalog=self._table)[0]
        names = self._set_names(table)

        table.add_column(Column(data=names, name='ID'))
        table.add_column(Column(data=['mag']*len(table), name='unit'))
        table.add_column(Column(data=['2013AJ....145...44Z']*len(table), name='bib'))

        return table

class DENIS_Catalog(Catalog):
    def __init__(self):
        Catalog.__init__(self, 'DENIS', vizier_table='B/denis', vizier_query=True)
        Vizier.ROW_LIMIT = -1

    def _set_names(self, table, **kwargs):
        return  ["DENIS %s" % i['DENIS'] for i in table]

    def _keys(self, filter=None):
        if filter is None:
            log.warn("No filter specified, loading just astrometric data.")
            return {'id_key':'ID',
                    'ra_key':'RAJ2000',
                    'dec_key':'DEJ2000'}

        if filter not in {'I', 'J', 'K'}:
            raise ValueError("DENIS catalog only have the following filters: I J K")

        return {'id_key':'ID',
                'ra_key':'RAJ2000',
                'dec_key':'DEJ2000',
                'flux_key':'%smag' % filter,
                'flux_error_key':'e_%smag' % filter,
                'flux_bib_key':'bib',
                'flux_unit_key':'unit'}

    def _load(self, center, radius, filter=None):
        table = Vizier.query_region(center, radius=Angle(radius), catalog=self._table)[0]
        names = self._set_names(table)

        table.add_column(Column(data=names, name='ID'))
        table.add_column(Column(data=['mag']*len(table), name='unit'))
        table.add_column(Column(data=['2005yCat.2263....0T']*len(table), name='bib'))

        return table

################################################################################
class CatalogLoader(object):
    '''
    This class helps to load a big list of catalogs to photobject colection.
    '''
    def __init__(self, *args, **kwargs):
        self._catalog_list = {'Simbad':Simbad_Catalog(),
                              'UCAC4':UCAC4_Catalog(),
                              'DENIS':DENIS_Catalog()}

    def query_catalog(self, catalog_name, center, radius, filter=None, *args, **kwargs):
        '''
        Load objects from a catalog.
        '''
        if catalog_name in self._catalog_list.keys():
            cat = self._catalog_list[catalog_name]
        else:
            raise ValueError('Catalog named %s not found. Use the list_catalogs() method tu see the available catalogs.' % catalog_name)

        return QueryData(cat.query(center, radius, filter, **kwargs), cat.get_keys(filter))

    def list_catalogs(self):
        '''
        List the available catalogs.
        '''
        for i in self._catalog_list.keys():
            print(i + ' : ' + self._catalog_list[i].description)
