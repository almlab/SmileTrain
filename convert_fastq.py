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
import util

def convert_entry(entry, check_id_parse=True):
    at_line, seq_line, qua_line = entry
    
    if check_id_parse:
        # try parsing the id
        rid = util.fastq_at_line_to_id(at_line)
        
    return [at_line, seq_line, util.illumina13_quality_to_18(qua_line)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a fastq from Illumina 1.3-1.7 format pipeline format')
    parser.add_argument('input', help='input fastq')
    parser.add_argument('output', help='output fastq')
    args = parser.parse_args()
    
    with open(args.input) as i:
        with open(args.output, 'w') as o:
            for entry in util.fastq_iterator(i):
                o.write(util.fastq_entry_list_to_string(convert_entry(entry)) + "\n")