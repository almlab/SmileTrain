'''
Script that removes primers from reads in fastq files.
'''

import argparse, itertools, sys
import usearch_python.primer_tools
from util import *

def mismatches(seq, primer, w):
    '''
    Calculate mismatches between a sequence and primer with window size w
    '''

    I = 0
    D = len(seq)
    for i in range(w):
        d = usearch_python.primer_tools.MatchPrefix(seq[i:], primer)
        if d < D:
            I = i
            D = d
    return [I, D]

def remove_primers(fastq, primer, max_primer_diffs):
    '''
    Remove well-matched primers from sequences in a fastq file

    Parameters
    fastq : sequence or iterator of strings
        lines in the fastq file
    primer : string
        the primer sequence to be removed
    max_primer_diffs : int
        the maximum number of mismatches allowed before throwing out the read

    Returns nothing. All output is delivered with 'print'.
    '''

    primer_length = len(primer)

    # loop over sets of 4 lines at a time
    for at_line, seq_line, plus_line, quality_line in itertools.izip(*[iter(fastq)] * 4):
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

        # find the primer position in the sequence
        primer_start_index, primer_diffs = mismatches(seq_line, primer, 15)

        # if we don't find a good match, move on. otherwise, trim the sequence and the
        # quality line and print a new fastq entry.
        if primer_diffs <= max_primer_diffs:
            primer_end_index = primer_start_index + primer_length
            trimmed_seq = seq_line[primer_end_index:]
            trimmed_quality = quality_line[primer_end_index:]
            print "\n".join([at_line, trimmed_seq, '+', trimmed_quality])


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help='input fastq file')
    parser.add_argument('primer', help='primer sequence')
    parser.add_argument('-m', '--max_primer_diffs', default=0, type=int, help='maximum number of nucleotide mismatches in the primer (default: 0)')
    args = parser.parse_args()

    remove_primers(open(args.fastq), args.primer, args.max_primer_diffs)