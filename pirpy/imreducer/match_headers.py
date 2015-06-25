#The list of the header keys to check if the instrument configure is the same in 2 images
opd_checkkeys = ['INSTRUME','PREAMP','HBIN','VBIN','EXPTIME','FILTER',
                 'GAIN','VSHIFT','IMGRECT','SUBRECT']

def match_header_fields(header1, header2, keys):
    '''
    This routine check if some fields of 2 fits headers matches.

    Parameters:
        header1 : dict
            Header to check
        header2 : dict
            Header to check
        keys : list of strings
            List of the keys to compare in headers

    Returns:
        bool: True if the selected valeus matches, False otherwise.
    '''
    for i in keys:
        if header1[i] != header2[i]:
            print('Key '+ i + ' don\'t match: ' +
                  str(header1[i]) + '     ' + str(header2[i]))
            return False
    return True

