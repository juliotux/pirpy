'''
Implements the photometry storage for a single object.
'''

from __future__ import division

from astropy.table import Table
from astropy.io import ascii, fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import Angle
from astropy.stats import median_absolute_deviation as mad
import numpy as np
from ast import literal_eval

from ..math.list_tools import to_list, match_lengths
from ..math.mag_tools import *
from ..math.stats import nanmad, nanmad_std
from .object_catalog import *
from ..mp import mult_ret

from ..log import log

__all__ = ['ResultStore', 'Frame']

phot_type = np.dtype([('id', 'a32'), ('x', 'float32'), ('y', 'float32'), ('flux', 'float64'), ('flux_err', 'float64')])
lc_type = np.dtype([('jd', 'float64'), ('relative_flux', 'float64')])

def params_to_sting(params):
    final = np.zeros(len(params), dtype='a256')
    for i in range(len(params)):
        final[i] = str(params[i])
    return final

class Frame(object):
    '''
    Stores the result of a single frame photometry.
    '''
    def __init__(self, jd):
        self.jd = jd
        self.data = np.zeros(0, dtype=phot_type)
        self.params = None
        self.params_errors = None

    def add_result(self, id, x, y, flux, flux_err, params=None, params_errors=None):
        self.add_results([id],[x],[y],[flux],[flux_err],[params],[params_erros])

    def add_results(self, id, x, y, flux, flux_err, params=None, params_errors=None):
        a = np.array(zip(id, x, y, flux, flux_err), dtype=phot_type)
        self.data = np.append(self.data, a)
        if params is not None:
            if self.params is None:
                self.params = np.zeros(0, dtype=params.dtype)
            self.params = np.append(self.params, params)
        if params_errors is not None:
            if self.params_errors is None:
                self.params_errors = np.zeros(0, dtype=params.dtype)
            self.params_errors = np.append(self.params_errors, params_errors)

    def rename_object(self, old_id, new_id):
        for i in range(len(self.data)):
            if self.data[i]['id'] == old_id:
                self.data[i]['id'] = new_id

    def compare(self, id1, id2):
        return self.data[self.data['id'] == id1]['flux']/self.data[self.data['id'] == id2]['flux']

    def instrumental_magnitude(self, id):
        return -2.5*np.log10(self.data[self.data['id'] == id]['flux'])

    def instrumental_magnitudes(self):
        return self.data['id'], -2.5*np.log10(self.data['flux'])

class ResultStore(object):
    '''
    Stores the result of various frames (possible grouped).
    '''
    def __init__(self, filter):
        self.frames = {} #Store the frames with a name and a Frame class
        self.groups = {} #Store the arrays containing the frame names of the groups

        self.filter = filter
        self.object_catalog = ObjectCatalog(self.filter)

    def _objects_in_frame(self, frame_name):
        return self.frames[frame_name].data['id']

    def _add_frame_to_groups(self, frame, groups=None):
        groups=to_list(groups)

        for i in groups:
            if i is not None:
                if i in self.groups.keys():
                    self.groups[i] = self.groups[i].union(set([frame]))
                else:
                    self.groups[i] = set([frame])

    def _get_frame_groups(self, frame):
        groups = []
        for i in self.groups.keys():
            if frame in self.groups[i]:
                groups.append(i)
        return groups

    def rename_star(self, old_name, new_name):
        '''
        Renames a star in the list.
        '''
        if new_name is None:
            log.warn('Object %s will not be renamed to None.' % old_name)
        else:
            if not old_name in self.object_catalog.get_ids():
                raise ValueError('Star %s not in the database' % old_name)
            if new_name != old_name:
                self.object_catalog.rename_star(old_name, new_name)
                for i in self.frames.keys():
                    self.frames[i].rename_object(old_name, new_name)

    def match_point(self, ra, dec):
        return self.object_catalog.match_point(ra, dec)

    def add_results(self, frame, jd, id, x, y, flux, flux_err, params=None, params_errors=None, ra=None, dec=None, groups=None):
        self.frames[frame] = Frame(jd)
        id = to_list(id)
        x = to_list(x)
        y = to_list(y)
        flux = to_list(flux)
        flux_err = to_list(flux_err)

        id, x, y, flux, flux_err = match_lengths([id, x, y, flux, flux_err], len(id))

        self.frames[frame].add_results(id, x, y, flux, flux_err, params, params_errors)

        self._add_frame_to_groups(frame, groups)

    def light_curve(self, obj, ref, group=None):
        lc = np.zeros(0, dtype=lc_type)

        if group is None:
            frames = self.frames.keys()
        else:
            frames = self.groups[group]

        for i in frames:
            if obj in self._objects_in_frame(i) and ref in self._objects_in_frame(i):
                lc = np.append(lc, np.array((self.frames[i].jd, self.frames[i].compare(obj, ref)), dtype=lc_type))

        return lc

    def diferential_magnitude(self, id, frame_group=None, n_it=100, mag_limits=(-100,100),
                               algorithm='montecarlo', percentage=0.5,
                               return_array=False):
        '''
        Calculate the median magnitude of an object, relative to all the others.]

        algorithm : 'montecarlo', 'median'
        '''
        frames = []
        if frame_group is None:
            for i in self.frames.keys():
                if id in self._objects_in_frame(i).flat:
                    frames.append(i)
        else:
            for i in self.groups[frame_group]:
                if id in self._objects_in_frame(i).flat:
                    frames.append(i)

        if algorithm == 'montecarlo':
            obj_list = set([])
            for i in frames:
                obj_list = obj_list.union(self._objects_in_frame(i).flat)
            try:
                obj_list.remove(id)
            except:
                dummy = 0

            objs = np.empty(0, dtype='a32')
            for i in obj_list:
                mag = self.object_catalog.get_mag(i)
                if not np.isnan(mag) and mag_limits[0] <= mag <= mag_limits[1]:
                    objs = np.append(objs, np.array([i], dtype='a32'))

            del obj_list

            results = np.array([np.nan]*n_it, dtype='float64')
            n_it_objs = int(len(objs)*percentage)

            for i in range(n_it):
                it_objs = np.random.choice(objs, n_it_objs)

                it_mags = np.zeros(0, dtype='float64')

                for k in frames:
                    f = self.frames[k]
                    f_mags = np.array([np.nan]*n_it_objs, dtype='float64')

                    obj_mag = f.instrumental_magnitude(id)
                    if len(obj_mag) > 0:
                        obj_mag = obj_mag[0]
                        for j in range(n_it_objs):
                            ref_mag = f.instrumental_magnitude(it_objs[j])
                            if len(ref_mag) > 0:
                                ref_mag = ref_mag[0]
                                ref_cat_mag = self.object_catalog.get_mag(it_objs[j])
                                f_mags[j] = obj_mag + (ref_cat_mag - ref_mag)
                            else:
                                f_mags[j] = np.nan

                    it_mags = np.append(it_mags, np.nanmedian(f_mags))
                results[i] = np.nanmedian(it_mags)

        elif algorithm == 'median':
            results = np.array([np.nan]*len(frames), dtype='float64')

            for i in range(len(frames)):
                mag_id, mags = self.frames[frames[i]].instrumental_magnitudes()
                obj_mag = mags[mag_id == id][0]
                cat_mags = np.zeros_like(mags)
                for j in range(len(mags)):
                    cat_mags[j] = self.object_catalog.get_mag(mag_id[j])

                results[i] = np.nanmedian(obj_mag + (cat_mags - mags))

        if return_array:
            return results
        else:
            return np.nanmedian(results), nanmad_std(results)

    def save_fits(self, fname, group=None):
        '''
        Saves the photometry frames in a fits file, including the object catalog.
        '''
        if group is None:
            frames = self.frames.keys()
        else:
            frames = self.groups[group]

        #creates the hdus for all the frames and index it
        fi = np.dtype([('index','int16'),('name','S64')])
        hdus = [None, None, None]
        frames_id = np.zeros(0, dtype=fi)
        frame_id = 2
        for i in frames:
            frame_id = frame_id + 1
            frames_id = np.append(frames_id, np.array((frame_id, i), dtype=fi))
            hdr = fits.Header()
            hdr['NAME'] = i
            hdr['JD'] = self.frames[i].jd
            for j,k in zip(self._get_frame_groups(i),range(len(self._get_frame_groups(i)))):
                hdr['GROUP%i' % k] = j
            columns = []
            columns.append(fits.Column(name='id', format='32A', array=self.frames[i].data['id']))
            columns.append(fits.Column(name='x', format='D', array=self.frames[i].data['x']))
            columns.append(fits.Column(name='y', format='D', array=self.frames[i].data['y']),)
            columns.append(fits.Column(name='flux', format='E', array=self.frames[i].data['flux']),)
            columns.append(fits.Column(name='flux_error', format='E', array=self.frames[i].data['flux_err']))

            if self.frames[i].params is not None:
                hdr['PARAMS'] = str(self.frames[i].params.dtype)
                columns.append(fits.Column(name='fit_parameters', format='256A', array=params_to_sting(self.frames[i].params)))
            if self.frames[i].params_errors is not None:
                hdr['PARAMS'] = str(self.frames[i].params.dtype)
                columns.append(fits.Column(name='fit_parameters_errors', format='256A', array=params_to_sting(self.frames[i].params_errors)))
            hdus.append(fits.BinTableHDU.from_columns(columns, header=hdr))

        hdus[2] = fits.BinTableHDU.from_columns(
                        [fits.Column(name='id', format='32A', array=self.object_catalog._cat['id']),
                         fits.Column(name='ra', format='D', array=self.object_catalog._cat['ra']),
                         fits.Column(name='dec', format='D', array=self.object_catalog._cat['dec']),
                         fits.Column(name='mag', format='E', array=self.object_catalog._cat['mag']),
                         fits.Column(name='mag_error', format='E', array=self.object_catalog._cat['mag_err']),
                         fits.Column(name='mag_unit', format='8A', array=self.object_catalog._cat['mag_unit'])])

        hdus[1] = fits.BinTableHDU.from_columns([fits.Column(name='index', format='I', array=frames_id['index']),
                                                 fits.Column(name='name', format='64A', array=frames_id['name'])])
        hdus[0] = fits.PrimaryHDU()

        hdulist = fits.HDUList(hdus)
        hdulist.writeto(fname)

    def load_fits(self, fname):
        f = fits.open(fname)

        #loading the object catalog
        self.object_catalog._cat = np.zeros(0, dtype=coord_type)
        self.object_catalog._cat = np.append(self.object_catalog._cat, np.array(zip(
            f[2].data['id'], f[2].data['ra'], f[2].data['dec'], f[2].data['mag'],
            f[2].data['mag_error'], f[2].data['mag_unit']
        ), dtype=coord_type))

        #load the frames
        for i,name in zip(f[1].data['index'], f[1].data['name']):
            fd = f[i].data
            fh = f[i].header

            if 'fit_parameters' in fd.columns.names:
                fit_parameters = [None]*len(fd['fit_parameters'])
                for j in range(len(fd['fit_parameters'])):
                    try:
                        fit_parameters[j] = literal_eval(fd['fit_parameters'][j])
                    except ValueError:
                        fit_parameters[j] = literal_eval(fd['fit_parameters'][j].replace('nan','\'nan\'').replace('inf','\'inf\''))
                fit_parameters = np.array(fit_parameters, dtype=literal_eval(fh['PARAMS']))
            else:
                fit_parameters = None

            if 'fit_parameters_errors' in fd.columns.names:
                fit_parameters_errors = [None]*len(fd['fit_parameters_errors'])
                for j in range(len(fd['fit_parameters_errors'])):
                    try:
                        fit_parameters_errors[j] = literal_eval(fd['fit_parameters_errors'][j])
                    except ValueError:
                        fit_parameters_errors[j] = literal_eval(fd['fit_parameters_errors'][j].replace('nan','\'nan\'').replace('inf','\'inf\''))
                fit_parameters_errors = np.array(fit_parameters_errors, dtype=literal_eval(fh['PARAMS']))
            else:
                fit_parameters_errors = None

            self.add_results(name, fh['JD'], fd['id'], fd['x'], fd['y'], fd['flux'], fd['flux_error'],
                             params=fit_parameters, params_errors=fit_parameters_errors)

            for k in fh:
                if 'GROUP' in k:
                    self._add_frame_to_groups(name, fh[k])
