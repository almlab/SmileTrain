#!/usr/bin/env python

'''
Convert qiime labels to smile train labels

e.g., ">samplename_123 lots of garbage stuff" -> ">sample=name;123"
'''

import argparse, re, sys
from Bio import SeqIO
from SmileTrain import util

def convert_qiime_label(record):
    ''' ">samplename_123 lots of garbage stuff" -> ">samplename;123" '''
    # split into "samplename" and "123"
    m = re.match("(.+)_(\d+)", record.id)
    if m is None: raise RuntimeError("record id '%s' did not parse" %(sid))
    sample_name, sequence_number = m.groups()
        
    # create a new id and throw out the other junk
    record.id = "sample=%s;%s" %(sample_name, sequence_number)
    record.description = ''
        
    return record
    

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output fasta (default stdout)')
    args = parser.parse_args()
    
    for record in SeqIO.parse(args.fasta, 'fasta'):
        SeqIO.write(convert_qiime_label(record), args.output, 'fasta')