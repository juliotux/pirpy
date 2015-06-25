'''
Functions to easy generate lists of fits files names, following some standarts.
'''

from os.path import join

def from_basename(basename, start_count, n_img, pad, sep='_', base_dir=''):
    '''
    This routine generates a file list, based on a basename, incrementing from
    a start_count by n_img times. The names follows the standart of the
    aquisition program from Pico dos Dias Observatory, and looks like:

        basename_0000.fits
        basename_0001.fits
        basename_0002.fits

    As this kind of name is used for other observatories, with different
    basename-number separator, you can set the sep variable to customize it.

    Args:
        basename : string
            The basename of the file list.
        start_count : int
            The number of the first image.
        n_img : int
            The total number of images.
        pad : int
            The number of digits to increment, like 4 for 0000.fits
        sep (optional) : string
            The string that separates the basename and the number. Default: '_'
        base_dir (optional) : string
            The path to be appended before the name.

    Returns:
        list of string: A list with the filenames.
    '''
    return [join(base_dir,basename + sep + str(i+start_count).rjust(pad,'0') + '.fits')
            for i in xrange(n_img)]

def from_file(filename):
    '''
    Generates a list of files from a file that contains the list. What???
    Example:

        $ ls *.fits >> mylist.txt

        mylist.txt
        ----
        file01.fits
        file02.fits
        file03.fits
        ----

    This code will read this and return:

        ['file01.fits','file02.fits','file03.fits']

    Parameters:
        filename : string
            The file that contains the list.

    Returns:
        list of string: The list of the filenames found inside the file.
    '''
    f1 = open(filename, 'rb').readlines()
    f2 = []
    for i in f1:
        f2.append(i.strip('\n'))
    return f2
