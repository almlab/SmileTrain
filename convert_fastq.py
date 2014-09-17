#!/usr/bin/env python

'''
Converts a fastq file from Illumina 1.3-1.7 format to a format that will work with this
pipeline, specifically:
    * an @ line that matches 1.4-1.7 format @stuff:with:colons#BARCODE/1 (or /2)
    * quality line using Illumina 1.8 (base 33) quality scores
'''

import argparse, sys
from Bio import SeqIO


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='convert a fastq from Illumina 1.3-1.7 format to 1.8 format', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fastq', help='input')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='output fastq (default: stdout)')
    args = parser.parse_args()
    
    for record in SeqIO.parse(args.fastq, 'fastq-illumina'):
        record.description = ''
        SeqIO.write(record, args.output, 'fastq')