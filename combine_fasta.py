#!/usr/bin/env python

'''
Combine multiple fasta files.
'''

import argparse
from Bio import SeqIO
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine multiple fasta files')
    parser.add_argument('fasta', nargs='+', help='input fasta')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    args = parser.parse_args()

    for fasta in args.fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            SeqIO.write(record, args.output, 'fasta')