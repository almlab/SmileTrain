#!/usr/bin/env python

'''
Get the entries from the dereplicated fasta file that have sequences that did not hit the
reference-based search.
'''

import uc2otus, util, argparse, sys, re

def missed_labels(uc_lines):
    '''uc lines -> labels for missed sequences'''
    
    labels = []
    for line in uc_lines:
        hit, label, otu = uc2otus.parse_uc_line(line)
        if hit == 'N': labels.append(label)
            
    return labels

def matching_fasta_entries(labels, fastas):
    new_fasta = ""
    for label, seq in fastas:
        if label in labels:
            new_fasta += ">%s\n%s\n" % (label, seq)
    
    return new_fasta
            
            

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('uc', help='input uc file')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output fasta file (default stdout)')
    args = parser.parse_args()
    
    with open(args.uc) as f:
        labels = missed_labels(f)
        
    with open(args.fasta) as f:
        fastas = util.fasta_entries(f, output_type='list')
        new_fasta = matching_fasta_entries(labels, fastas)
    
    # rename size in label
    new_fasta = re.sub('counts', 'size', new_fasta)
    args.output.write(new_fasta)
