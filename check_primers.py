#!/usr/bin/env python

'''
Script that checks for the presence of a primer in a fasta.
'''

import argparse, itertools, sys, tempfile
import util, util_primer


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help='input fastq file')
    parser.add_argument('primer', help='primer sequence')
    parser.add_argument('-l', '--log_frac', type=int, help='negative logarithm of fraction of sequences to take (default: 4')
    parser.add_argument('-m', '--max_primer_diffs', default=0, type=int, help='maximum number of nucleotide mismatches in the primer (default: 0)')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    
    args = parser.parse_args()

    skip = int(10 ** -args.log_frac)
    r = util_primer.PrimerRemover(open(args.fastq), args.primer, args.max_primer_diffs, skip=skip)
    r.check_entries()
    
    args.output.write(r.diagnostic_message())
