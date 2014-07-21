#!/usr/bin/env python

'''
Check if a primer (or its reverse complement) actually appears in the target fastq file
'''

import argparse, re, sys
from SmileTrain import util, util_primer
import Bio.Seq


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input fastq')
    parser.add_argument('primer', help='input primer')
    parser.add_argument('-c', '--reverse_complement', action='store_true', help='check for reverse complement instead?')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output fasta (default stdout)')
    args = parser.parse_args()
    
    with open(args.input) as f:
        if args.reverse_complement:
            primer = Bio.Seq.Seq(args.primer).reverse_complement().tostring()
        else:
            primer = args.primer
            
        p = util_primer.PrimerRemover(f, primer, 0)
        p.check_entries()
        print p.diagnostic_message()