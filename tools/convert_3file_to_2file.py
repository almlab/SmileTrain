#!/usr/bin/env python

'''
Convert the three-file raw Illumina data (forward, reverse, and index reads) into the
two-file format used by SmileTrain.
'''

import argparse, sys, os
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from SmileTrain import util

def lid(fastq_entry):
    '''get the rid and convert "@M01056:blabla:1382 2:N:0:0" -> "@M01056:blabla:1382"'''
    return fastq_entry[0].split()[0]

def new_id_line(location, barcode, direction):
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
    
    with open(args.forward_in) as fi, open(args.reverse_in) as ri, open(args.index) as index, open(args.forward_out, 'w') as fo, open(args.reverse_out, 'w') as ro:
        fi_fq = util.fastq_iterator(fi)
        ri_fq = util.fastq_iterator(ri)
        i_fq = util.fastq_iterator(index)
        
        for fi_entry in fi_fq:
            ri_entry = ri_fq.next()
            i_entry = i_fq.next()
            
            # check that we are looking at the same location
            assert(lid(fi_entry) == lid(ri_entry) == lid(i_entry))
            
            # get new location id
            location = lid(fi_entry)
            barcode = i_entry[1]
            
            # make new forward and reverse entries
            fo_entry = [new_id_line(location, barcode, "1"), fi_entry[1], fi_entry[2]]
            ro_entry = [new_id_line(location, barcode, "2"), ri_entry[1], ri_entry[2]]
            
            # output the entries to their new places
            fo.write(util.fastq_entry_list_to_string(fo_entry) + "\n")
            ro.write(util.fastq_entry_list_to_string(ro_entry) + "\n")