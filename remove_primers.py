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

If a reverse primer is specified, that sequence is removed from the end of the sequence.
Because merging causes the end of the sequence to be the reverse complement of the
reverse primer, it is actually the rev comp of the rev primer that will be removed.

Entries that do not match the primer(s) are dropped.
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
    parser.add_argument('-q', '--reverse_primer', default=None, help='reverse primer sequence')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output fastq (default stdout)')
    args = parser.parse_args()

    r = util_primer.PrimerRemover(args.fastq, args.primer, args.max_primer_diffs, reverse_primer=args.reverse_primer, out=args.output)
    r.print_entries()

    if args.log is not None:
        with open(args.log, 'w') as f:
            f.write(r.diagnostic_message())