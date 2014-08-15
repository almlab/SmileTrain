#!/usr/bin/env python

'''
Converts a fastq file from Illumina 1.3-1.7 format to a format that will work with this
pipeline, specifically:
    * an unchanged @ line that matches 1.4-1.7 format @stuff:with:colons#BARCODE/1 (or /2)
    * sequence line
    * a line with just +
    * quality line using Illumina 1.8 (base 33) quality scores
'''

import itertools, os, os.path, sys, argparse, itertools, shutil
from Bio import SeqIO
import util

def convert_record(record, offset=-31):
    '''
    Convert a BioPython fastq record
    
    record : Seq
        to be converted
    offset : int
        how to change the quality score (default: -31, which goes from Illumina 1.3 to 1.8)
        
    returns : Seq
    '''
    
    scores = record.letter_annotations['phred_quality']
    new_scores = [s + offset for s in scores]
    record.scores = new_scores
    return record

def convert_record_illumina13_to_18(record):
    return convert_record(record, offset=-31)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a fastq from Illumina 1.3-1.7 format to 1.8 format')
    parser.add_argument('input', help='input fastq')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='output fastq (default: stdout)')
    args = parser.parse_args()
    
    for record in SeqIO.parse(args.input, 'fastq'):
        record = convert_record_illumina13_to_18(record)
        SeqIO.write(record, args.output, 'fastq')