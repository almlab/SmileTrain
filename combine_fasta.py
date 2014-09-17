#!/usr/bin/env python

'''
Combine multiple fasta files.
'''

import argparse, sys
from Bio import SeqIO
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine multiple fasta files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', nargs='+', help='input fasta')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='output file')
    args = parser.parse_args()

    for fasta in args.fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            SeqIO.write(record, args.output, 'fasta')