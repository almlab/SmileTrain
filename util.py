import re, string, sys, time, itertools, os
import usearch_python.primer

def fasta_entries(lines):
    '''
    Yield [id, sequence] pairs from a fasta file. Sequence allowed to run over multiple lines.

    lines : sequence or iterator of strings
        lines from the fasta file

    Yields [id (without the >), sequence] pairs
    '''

    sid = None
    sequence = ''
    for line in lines:
        line = line.rstrip()

        # check if this is the first line
        if sid is None:
            assert(line.startswith('>'))
            sid = line[1:]
        else:
            if line.startswith('>'):
                assert(sid != '')
                assert(sequence != '')

                yield [sid, sequence]
                sid = line[1:]
                sequence = ''
            else:
                sequence += line

    yield [sid, sequence]

def fastq_iterator(lines, check_sigils=True, check_lengths=True, output_type='list'):
    '''
    Yield [at line, seq line, quality line] entries from a fastq file. All lines are right-
    trimmed, and the plus line (line 3) is dropped.

    lines : sequence or iterator of strings
        lines from the fastq file
    check_sigils : bool (default true)
        assert that the first line starts with @ and the third with +?
    check_lengths : bool (default true)
        assert that the sequence and quality lines are the same length
    output_type : string 'string' or 'list' (default 'list')
        output format

    Returns : if 'string', return a single string; if 'list', return a list [at line, sequence line, quality line]
    '''

    for at_line, seq_line, plus_line, quality_line in itertools.izip(*[iter(lines)] * 4):
        # chomp all newlines
        at_line = at_line.rstrip()
        seq_line = seq_line.rstrip()
        plus_line = plus_line.rstrip()
        quality_line = quality_line.rstrip()

        # check that the two lines with identifiers match our expectations
        if check_sigils:
            assert(at_line.startswith('@'))
            assert(plus_line.startswith('+'))

        # check that the sequence and quality lines have the same number of nucleotides
        if check_lengths:
            assert(len(seq_line) == len(quality_line))

        entry = [at_line, seq_line, quality_line]
        if output_type == 'list':
            yield entry
        elif output_type == 'string':
            yield fastq_entry_list_to_string(entry)

def fastq_entry_list_to_string(entry):
    '''[at line, sequence line, quality line] -> new-line separated fastq entry'''
    at_line, seq_line, quality_line = entry
    return "\n".join([at_line, seq_line, '+', quality_line])

def fastq_at_line_to_id(line):
    '''"@lolapolooza/1" -> "lolapolooza"'''
    m = re.match('@(.+)/[12]', line)
    if m is None:
        raise RuntimeError("fastq at line did not parse: %s" % line)
    else:
        rid = m.group(1)
        
    return rid

# Note that the 1.8 code J (score 41) has not equivalent in 1.3, so I map it to h (score 40)
illumina13_codes = '''ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghh'''
illumina18_codes = '''"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ'''
illumina_code_table = string.maketrans(illumina13_codes, illumina18_codes)

def illumina13_quality_to_18(quality_line):
    '''convert Illumina 1.3-1.7 quality scores to Illumina 1.8 quality scores'''
    return quality_line.translate(illumina_code_table)

def mismatches(seq, primer, w):
    '''
    Calculate mismatches between a sequence and primer with window size w.
    Returns the starting index and number of mismatches for the best match.
    '''

    I = 0
    D = len(seq)
    for i in range(w):
        d = usearch_python.primer.MatchPrefix(seq[i:], primer)
        if d < D:
            I = i
            D = d
    return [I, D]

def listify(inp):
    '''Ensure that input is a list
    
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
        tests = [os.path.isfile(filename) for filename in filenames]
        if False in tests:
            bad_names = " ".join([filename for filename, test in zip(filenames, tests) if test == False])
            raise RuntimeError("input file(s) missing: %s" % bad_names)
    
def check_for_nonempty(filenames, dry_run=False):
    '''assert that each file exists and is nonempty'''
    
    filenames = listify(filenames)
    
    if dry_run:
        print "dry run: test that files are non-empty: " + " ".join(filenames)
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
        print "dry run: test that destinations are free: " + " ".join(filenames)
    else:
        tests = [os.path.isfile(filename) for filename in filenames]
        if True in tests:
            bad_names = " ".join([filename for filename, test in zip(filenames, tests) if test == True])
            raise RuntimeError("output file(s) already exist: %s" % bad_names)
    
def is_executable(filename):
    '''check if a filename exists and is executable'''
    return os.path.isfile(filename) and os.access(filename, os.X_OK)

def message(text, indent=2):
    # print message to stderr
    space = ' ' * indent
    text = re.sub('\n', '\n%s' %(space), text)
    sys.stderr.write('%s%s\n' %(space, text))


def error(text, indent=2):
    # print message to stderr and quit
    space = ' ' * indent
    text = re.sub('\n', '\n%s' %(space), text)
    sys.stderr.write('%s%s\n' %(space, text))
    quit()


def read_list(fn, dtype=str):
    # read file as list
    x = [dtype(line.rstrip()) for line in open(fn)]
    return x


def read_dataframe(fn, index_dtype=str, columns_dtype=str):
    import pandas as pd
    # read file as pandas dataframe
    x = pd.read_table(fn, sep='\t', header=0, index_col=0)
    x.index = x.index.astype(index_dtype)
    x.columns = x.columns.astype(columns_dtype)
    return x


def read_tseries(fn):
    import pandas as pd
    # read file as pandas time series
    return read_dataframe(fn, index_dtype=float, columns_dtype=str)


def iter_fst(fn):
    # generator that iterates through [sid, seq] pairs in a fasta file
    sid = ''
    seq = ''
    for line in open(fn):
        line = line.rstrip()
        if line.startswith('>'):
            if seq != '':
                yield [sid, seq]
            sid = line[1:]
            seq = ''
        else:
            seq += line
    yield [sid, seq]


def iter_fsq(fn):
    # generator that iterates through records in a fastq file
    record = []
    i = 0
    for line in open(fn):
        i += 1
        if i % 4 == 1:
            if len(record) > 0:
                yield record
            record = []
        record.append(line.rstrip())
    yield record


def read_fst(fn, reverse=False):
    # read fasta file as dictionary
    fst = {}
    for [sid, seq] in iterfst(fst):
        if reverse == False:
            fst[sid] = seq
        elif reverse == True:
            fst[seq] = sid
    return fst


def cycle(x):
    # an efficient way to cycle through a list (similar to itertools.cycle)
    while True:
        for xi in x:
            yield xi


class timer():
    # generator that measures elapsed time
    def __init__(self):
        self.t = [time.time()]
    def __iter__(self):
        return self
    def next(self):
        self.t.append(time.time())
        return self.t[-1] - self.t.pop(0)

rctab = string.maketrans('ACGTacgt','TGCAtgca')
def reverse_complement(x):
    return x[::-1].translate(rctab)


