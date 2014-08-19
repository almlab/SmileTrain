import re, string, sys, time, itertools, os, subprocess
import usearch_python.primer

def listify(inp):
    '''
    Ensure that input is a list
    
    inp : string or list
    
    returns : list
        if input is a list, return that; otherwise return [string]
    '''
    
    if isinstance(inp, str):
        return [inp]
    elif isinstance(inp, list):
        return inp
    else:
        raise RuntimeError("don't know how to stringify input: " + inp)

def check_for_existence(filenames, dry_run=False):
    '''assert that each of filenames does exist'''

    filenames = listify(filenames)

    if dry_run:
        print "dry run: test for existence of files: " + " ".join(filenames)
    else:
        # running ls first seems to prevent spurious empties
        subprocess.check_output(['ls', '-lah'])
        
        tests = [os.path.isfile(filename) for filename in filenames]
        if False in tests:
            bad_names = " ".join([filename for filename, test in zip(filenames, tests) if test == False])
            raise RuntimeError("file(s) missing: %s" % bad_names)
    
def check_for_nonempty(filenames, dry_run=False):
    '''assert that each file exists and is nonempty'''
    
    filenames = listify(filenames)
    
    if dry_run:
        message("dry run: test that files are non-empty: " + " ".join(filenames), indent=4)
    else:
        check_for_existence(filenames)
    
        tests = [os.stat(fn).st_size > 0 for fn in filenames]
        if False in tests:
            bad_names = " ".join([fn for fn, t in zip(filenames, tests) if t == False])
            raise RuntimeError("file(s) empty: " + bad_names)

def check_for_collisions(filenames, dry_run=False):
    '''assert that each of filenames does not exist'''

    # correct a string into a list if needed
    filenames = listify(filenames)

    if dry_run:
        message("dry run: test that destinations are free: " + " ".join(filenames), indent=4)
    else:
        tests = [os.path.isfile(filename) for filename in filenames]
        if True in tests:
            bad_names = " ".join([filename for filename, test in zip(filenames, tests) if test == True])
            raise RuntimeError("output file(s) already exist: %s" % bad_names)
    
def is_executable(filename):
    '''check if a filename exists and is executable'''
    return os.path.isfile(filename) and os.access(filename, os.X_OK)

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


banner = '''
 ________            _____ ______      ________                _____          
 __  ___/_______ ___ ___(_)___  /_____ ___  __/______________ ____(_)_______  
 _____ \ __  __ `__ \__  / __  / _  _ \__  /   __  ___/_  __ `/__  / __  __ \ 
 ____/ / _  / / / / /_  /  _  /  /  __/_  /    _  /    / /_/ / _  /  _  / / / 
 /____/  /_/ /_/ /_/ /_/   /_/   \___/ /_/     /_/     \__,_/  /_/   /_/ /_/ 

'''