#!/usr/bin/env python

'''
Check which barcodes appear in the fasta, and hope that their distributions are as
expected. Allow 0 mismatches.
'''

import argparse, sys, os, csv
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from SmileTrain import map_barcodes, util
from Bio import SeqIO


def count_barcodes(fastq, barcode_map):
    '''count the number of appearances of each barcode in a fastq'''
    
    # initialize the count dictionary
    counts = {sample: 0 for sample in barcode_map.values()}
    counts['mapped'] = 0
    counts['total'] = 0
    
    for record in SeqIO.parse(fastq, 'fastq'):
        barcode_read, read_direction = map_barcodes.parse_barcode(record)
        if barcode_read in barcode_map:
            sample = barcode_map[barcode_read]
            counts[sample] += 1
            counts['mapped'] += 1
        
        counts['total'] += 1
    
    return counts

def write_counts_report(counts, f):
    '''log the counts in a filehandle f'''
    
    # sort the results by abundance
    sorted_counts = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    
    percent = lambda n: "%.2f%%" %(float(n) / counts['total'] * 100)
    
    w = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    w.writerow(['sample', '# reads', '% reads'])
        
    for sample, n in sorted_counts:
        w.writerow([sample, str(n), percent(n)])


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help='input fastq')
    parser.add_argument('barcodes', help='input barcode mapping (or list of barcode)')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output log (default stdout)')
    args = parser.parse_args()
    
    with open(args.barcodes) as f:
        barcode_map = map_barcodes.barcode_file_to_dictionary(f)
    
    counts = count_barcodes(args.fastq, barcode_map)
    
    write_counts_report(counts, args.output)
