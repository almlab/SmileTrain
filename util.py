import re, string, sys, time, itertools, os, subprocess
import usearch_python.primer

def listify(inp):
    '''
    Ensure that input is a list or tuple
    
    inp : list or single item
    
    returns : list
        if input is a list, return that; otherwise return [item]
    '''
    
    if isinstance(inp, list) or isinstance(inp, tuple):
        return inp
    else:
        return [inp]

def check_for_collisions(filenames):
    '''assert that each of filenames does not exist'''
    # correct a string into a list if needed
    filenames = listify(filenames)
    
    tests = [os.path.isfile(filename) for filename in filenames]
    if True in tests:
        bad_names = " ".join([filename for filename, test in zip(filenames, tests) if test == True])
        raise RuntimeError("output file(s) already exist: %s" % bad_names)

def message(text, indent=2):
    '''print message to stderr'''
    space = ' ' * indent
    text = re.sub('\n', '\n%s' %(space), text)
    sys.stderr.write('%s%s\n' %(space, text))

def error(text, indent=2):
    '''print message to stderr and quit'''
    space = ' ' * indent
    text = re.sub('\n', '\n%s' %(space), text)
    sys.stderr.write('%s%s\n' %(space, text))
    quit()


class timer():
    '''generator that measures elapsed time'''
    def __init__(self):
        self.t = [time.time()]
    def __iter__(self):
        return self
    def next(self):
        self.t.append(time.time())
        return self.t[-1] - self.t.pop(0)