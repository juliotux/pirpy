'''
Tool to make a path tree with easy error handle.
'''
import errno
import os

def mkdir_p(path):
    '''
    Creates a directory tree recursivelly, checking if it exists.

    Parameters:
        path : string
            The path to be created.
    '''
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
