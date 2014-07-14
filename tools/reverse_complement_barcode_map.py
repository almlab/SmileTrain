#!/usr/bin/env python

'''
Take the reverse complements of the barcodes to make a new mapping file.
'''

import argparse, sys, csv
import Bio.Seq


def reverse_complement_map(m):
    '''[['sample1', 'AAA'], ...] -> [['sample1', 'TTT'], ...]'''
    return [[sample, Bio.Seq.Seq(code).reverse_complement().tostring()] for sample, code in m]


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('barcodes', help='input barcode mapping')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output mapping file (default stdout)')
    args = parser.parse_args()
    
    with open(args.barcodes) as f:
        r = csv.reader(f, delimiter='\t')
        old_map = [l for l in r]
        assert(len(old_map[0]) == 2)
        
    new_map = reverse_complement_map(old_map)
    
    w = csv.writer(args.output, delimiter='\t', quoting=csv.QUOTE_NONE)
    w.writerows(new_map)