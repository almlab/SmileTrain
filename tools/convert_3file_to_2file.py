#!/usr/bin/env python

'''
Convert the three-file raw Illumina data (forward, reverse, and index reads) into the
two-file format used by SmileTrain.
'''

import argparse, sys, os, itertools
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from Bio import SeqIO

def new_id(location, barcode, direction):
    '''(location, barcode, direction) -> location#barcode/direction'''
    return "%s#%s/%s" %(location, barcode, direction)

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('forward_in', help='input forward reads fastq')
    parser.add_argument('reverse_in', help='input reverse reads fastq')
    parser.add_argument('index', help='input index reads fastq')
    parser.add_argument('forward_out', help='output forward reads fastq')
    parser.add_argument('reverse_out', help='output reverse reads fastq')
    
    args = parser.parse_args()
    
    with open(args.forward_out, 'w') as fo, open(args.reverse_out, 'w') as ro:
        fi = SeqIO.parse(args.forward_in, 'fastq')
        ri = SeqIO.parse(args.reverse_in, 'fastq')
        ii = SeqIO.parse(args.index, 'fastq')
        
        for fr, rr, ir in itertools.izip(fi, ri, ii):
            # check that the reads match
            if not (fr.id == rr.id == ir.id):
                raise RuntimeError("ids in corresponding entries did not match: %s %s %s" %(fr.id, rr.id, ir.id))
            
            # the barcode is the sequence from the index read
            barcode = str(ir.seq)
            
            # amend the ids on the forward and reverse entries
            fr.id = new_id(fr.id, barcode, '1')
            rr.id = new_id(rr.id, barcode, '2')
            fr.description = ''
            rr.description = ''
            
            # write the entries to their output files
            SeqIO.write(fr, fo, 'fastq')
            SeqIO.write(rr, ro, 'fastq')