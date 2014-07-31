#!/usr/bin/env python

'''
Script that removes primers from reads in fastq files.

With a primer
    AAAAA

entries like
    @lolapolooza
    TAAAACATCATCATCAT
    +whatever
    abcdefghijklmnopq

become entries like
    @lolapolooza
    CATCATCATCAT
    +
    fghijklmnopq
'''

import argparse, itertools, sys
import util, util_primer


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help='input fastq file')
    parser.add_argument('primer', help='primer sequence')
    parser.add_argument('-m', '--max_primer_diffs', default=0, type=int, help='maximum number of nucleotide mismatches in the primer (default: 0)')
    parser.add_argument('-l', '--log', default=None, type=str, help='log file for successes, failures, and time elapsed')
    args = parser.parse_args()

    r = util_primer.PrimerRemover(open(args.fastq), args.primer, args.max_primer_diffs)
    r.print_entries()

    if args.log is not None:
        with open(args.log, 'w') as f:
            f.write(r.diagnostic_message())
