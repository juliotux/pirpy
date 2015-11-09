from multiprocessing import Pool
from subprocess import call

try:
    import os
    num_threads = os.environ['OMP_NUM_THREADS']
except:
    num_threads = 1

def mult_ret(func, arg, nprocess):
    '''
    Function to return a list of values, based on a list of arguments
    passed to a single function. Can be multiprocess or not.
    '''
    if nprocess > 1:
        pool = Pool(nprocess)
        return pool.map(func,arg)
    else:
        ret = [None]*len(arg)
        for i in range(len(arg)):
            ret[i]= func(arg[i])
        return ret

def call_mp(commands, nprocess):
    '''
    Use the subprocess.call standart function inside a multiprocessing.Pool.
    It's useful to work with external programs that takes some time to process.

    Parameters:
        commands : list of list of strings
            Basically, is the list of commands to pass to the subprocess.call
            function. Each command is a list of the program itself and it's
            arguments.
        nprocess : int
            The number of process to generate the Pool.
    '''
    if nprocess > 1:
        pool = Pool(processes=nprocess)
        pool.map(call,commands)
    else:
        for i in commands:
            call(i)
