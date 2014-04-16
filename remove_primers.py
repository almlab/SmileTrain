import argparse, primer, sys
from util import *

def mismatches(seq, subseq, w):
    # calculate mismatches between seq and subseq with window size w
    I = 0
    D = len(seq)
    for i in range(w):
        d = primer.MatchPrefix(seq[i:], subseq)
        if d < D:
            I = i
            D = d
    return [I, D]


def parse_args():
    # parse command line arguments
    parser = argparse.Argument_Parser()
    parser.add_argument('-i', default = '', help = 'input fasta file')
    parser.add_argument('-p', default = '', help = 'primer sequence')
    args = parser.parse_args()


def remove_primers():
    # remove primers



# Read input arguments

fasta = sys.argv[1]
prime = sys.argv[2]
MAX_PRIMER_DIFFS = int(sys.argv[3])
PL = len(prime)

# Remove primers
k = 0
for line in open(fasta):
    line = line.rstrip()
    k += 1
    if k%4 == 1:
        r = line[1:]
    elif k%4 == 2:
        seq = line
        I, D = mismatches(seq, prime, 15)
        if D > MAX_PRIMER_DIFFS:
            r = ''
        else:
            seq = seq[I+PL:]
    elif r != '' and k%4 == 0:
        qua = line[I+PL:]
        print '@%s\n%s\n+\n%s' %(r, seq, qua)
