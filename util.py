import re, string, sys, time, itertools
import usearch_python.primer

def fasta_entries(lines):
    '''
    Yield [id, sequence] pairs from a fasta file. Sequence allowed to run over multiple
    lines.

    Parameters
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

def fastq_iterator(lines, check_sigils=True, check_lengths=True):
    '''
    Yield [at line, seq line, quality line] entries from a fastq file. All lines are right-
    trimmed, and the plus line (line 3) is dropped.

    Parameters
    lines : sequence or iterator of strings
        lines from the fastq file
    check_sigils : bool (default true)
        assert that the first line starts with @ and the third with +?
    check_lengths : bool (default true)
        assert that the sequence and quality lines are the same length
    '''

    for at_line, seq_line, plus_line, quality_line in itertools.izip(*[iter(lines)] * 4):
        # chomp all newlines
        at_line = at_line.rstrip()
        seq_line = seq_line.rstrip()
        plus_line = plus_line.rstrip()
        quality_line = quality_line.rstrip()

        # check that the two lines with identifiers match our expectations
        assert(at_line.startswith('@'))
        assert(plus_line.startswith('+'))

        # check that the sequence and quality lines have the same number of nucleotides
        assert(len(seq_line) == len(quality_line))

        yield [at_line, seq_line, quality_line]

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


