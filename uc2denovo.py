#!/usr/bin/env python

'''
Get the entries from the dereplicated fasta file that have sequences that did not hit the
reference-based search.
'''

import uc2otus, util, argparse, sys, re
from Bio import SeqIO

def missed_labels(uc_lines):
    '''uc lines -> labels for missed sequences'''
    
    labels = []
    for line in uc_lines:
        hit, label, otu = uc2otus.parse_uc_line(line)
        if hit == 'N': labels.append(label)
            
    return labels

def matching_fasta_entries(labels, fasta):
    records = [record for record in SeqIO.parse(fasta, 'fasta') if record.id in labels]
    
    # rename 'counts' to 'size' to match usearch
    for record in records:
        record.id = re.sub('counts', 'size', record.id)
    
    return records
            

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('uc', help='input uc file')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output fasta file (default stdout)')
    args = parser.parse_args()
    
    with open(args.uc) as f:
        labels = missed_labels(f)
    
    records = matchind_fasta_entries(labels, SeqIO.parse(fasta, 'fasta'))
    SeqIO.write(records, args.output, 'fasta')