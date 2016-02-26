'''
Routine to wrapper the astrometry.net program.

This program calibrates fits images astrometrically, comparing the detected
sources with catalogs of stars. As the result, the program returns a new fits
file with WCS calibration intormation.
'''

import tempfile
import shutil
import os
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

from ..io.io_tool import mkdir_p
from ..mp import call_mp

__all__ = ['WCSCalib']

class WCSCalib():
    '''
    Instance to calibrate fits images with WCS. This class is basically
    a wrapper to the Astrometry.net software.

    Usage:
        ```
        #initialize an instance
        mywcscalib = WCSCalib('my/calibrated/dir',nprocess=4)
        mywcscalib.queue_files(myfilelist)
        mywcscalib.queue_files('singlefile.fits')
        mywcscalib.run()
        ```

    Arguments:
        calib_dir : string
            The location to output the result files.
        nprocess (optional) : int
            The number of parallel processes to use. Default: 1
        platescale_uncertain (optional) : int
            The uncertain in platescale from header, in platescale units
            (like 0.1 for 10%). Default: 0.1
    '''
    #TODO: Finish the documentation of the arguments
    def __init__(self, calib_dir, **kwargs):
        self.calib_dir = calib_dir
        self.args = kwargs
        self.queue = []

        self.nprocess = kwargs.get('nprocess',1)
        self.platescale_uncertain = kwargs.get('platescale_uncertain',0.1)
        self.do_guess_scale = kwargs.get('do_guess_scale',True)
        self.do_guess_coords = kwargs.get('do_guess_coords',True)
        self.search_radius = kwargs.get('search_radius',1)
        self.downsample = kwargs.get('downsample',1)
        self.no_plot = kwargs.get('no_plot',True)
        self.depth = kwargs.get('depth',20)
        self.keep_tree = kwargs.get('keep_tree',True) #To avoid conflicts

        self.guess_coords = kwargs.get('guess_coords',None) #{'ra':hourangle, 'dec': deg}
        self.guess_plate_scale = kwargs.get('guess_plate_scale',None) #approximated plate scale

    def queue_files(self, fname):
        '''
        Queue a file or file list

        Arguments:
            fname : string or list of strings
                The name(s) of the fits file(s) to calibrate.
        '''
        if isinstance(fname, basestring):
            self.queue.append(self.gen_commands(fname))
        elif isinstance(fname, list):
            for i in fname:
                self.queue.append(self.gen_commands(i))
        else:
            print('Error: Not a string or a list of strings.')

    def run(self):
        '''Run the queued commands.'''
        # Prepare the call
        dest_dirs = []
        tmp_dirs = []
        comms = []
        for i in self.queue:
            dest_dirs.append(i['dest_dir'])
            mkdir_p(i['dest_dir'])
            tmp_dirs.append(i['tmp_dir'])
            comms.append(i['commands'])

        # Call the software, with multiprocess
        call_mp(comms,self.nprocess)

        # Remove temporary folders
        for i in tmp_dirs:
            shutil.rmtree(i)

    def gen_commands(self, fname):
        '''
        Generates a list of commands to call astrometry.net

        Arguments:
            fname : string
                The path to the the FITS file to be aclibrated

        Returns:
            dict{'commands','tmp_dir','calib_dir'}
        '''

        # Allow user to choose if he want to keep the original tree or not
        destfile = None
        if self.keep_tree:
            destfile = fname
        else:
            destfile = os.path.basename(fname)

        dic = {'commands':['solve-field'],
               'tmp_dir':tempfile.mkdtemp(),
               'dest_dir':os.path.join(self.calib_dir,os.path.dirname(destfile))}

        coords = fieldsize_high = fieldsize_low = None

        if self.do_guess_scale or self.do_guess_coords:
            # Open FITS and make a guess of the center of the file, and get
            # informations about the field size. This makes astrometry.net fast
            f = fits.open(fname)
            if self.do_guess_coords:
                if self.guess_coords is not None:
                    a = self.guess_coords
                else:
                    a = f[0].header
                coords = SkyCoord(str(a['ra']).replace(',','.') + ' ' +
                                  str(a['dec']).replace(',','.'),
                                  unit=(u.hourangle, u.deg))
                dic['commands'].append('--dec=' + str(coords.dec.degree))
                dic['commands'].append('--ra=' + str(coords.ra.degree))
                dic['commands'].append('--radius=' + str(self.search_radius))


            if self.do_guess_scale:
                if self.guess_plate_scale is not None:
                    ps = self.guess_plate_scale
                else:
                    ps = float(f[0].header['PLATESCL'].strip(' '))
                dic['commands'].append('--scale-units')
                dic['commands'].append('arcsecperpix')
                dic['commands'].append('--scale-high')
                dic['commands'].append(str(ps * (1 + self.platescale_uncertain)))
                dic['commands'].append('--scale-low')
                dic['commands'].append(str(ps * (1 - self.platescale_uncertain)))
            f.close()

        if self.no_plot:
            dic['commands'].append('--no-plot')

        if self.downsample > 1 and isinstance(self.downsample,int):
            dic['commands'].append('--downsample')
            dic['commands'].append(str(self.downsample))

        dic['commands'].append('--no-fits2fits')
        dic['commands'].append('--depth')
        dic['commands'].append(str(self.depth))
        dic['commands'].append('--new-fits')
        dic['commands'].append(os.path.join(self.calib_dir,destfile))
        dic['commands'].append('--dir')
        dic['commands'].append(dic['tmp_dir'])
        dic['commands'].append(fname)

        return dic
