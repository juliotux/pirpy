'''
The allowed algorithms to use in the photometry.
    'photutils_aperture' : photutils astopy filiated package with aperture
                           photometry.
    'photutils_psf' : photutils astopy filiated package with PSF/PRF photometry.
    'sextractor' : sextractor algorithm implemented via SEP (Source Extraction
                   and Photometry) package (aperture).
'''

allowed = set(['photutils_aperture',
               'sextractor'])

todo = set(['photutils_psf'])
