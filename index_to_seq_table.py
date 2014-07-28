#!/usr/bin/env python

'''
Create a sequence table from the index file. Option to use the dereplicated fasta
for ordering by abundance.
'''

import re, sys, argparse
import util, util_index

def trimmed_fasta_entries(lines):
    '''create a list [[sequence id (without counts), sequence], ...]'''
    return [[util_index.parse_seq_sid(sid), seq] for sid, seq in util.fasta_entries(lines)]

def sparse_seq_count_table(index_lines):
    '''
    Create a sparse sequence count table using sample-sequence-abundance index.
    
    index_lines : list or iterator of strings
        lines from the index file (tab-separated sample, sequence ID, abundance)
        
    returns : dictionary of dictionaries
        {sample => {sequence_id => abundance, ...}, ...}
    '''
    
    table = {}
    for line in index_lines:
        sample, seq, abund = util_index.parse_index_line(line)
        
        if sample not in table:
            table[sample] = {seq: abund}
        else:
            if seq not in table[sample]:
                table[sample][seq] = abund
            else:
                table[sample][seq] += abund
                
    return table

def seq_table(table, lab_seq, samples=None):
    '''
    Output lines of a seq table.
    
    table : dictionary of dictionaries
        {sample => {seq_label => abundance, ...}, ...}
    lab_seq : list of lists
        [[sequence label, sequence], ...]
    samples : list or iterator of strings (default None)
        sample names in column order; or just make a new sorted list
        
    yields : strings
        lines in the otu table
    '''
    
    # make our own samples list if necessary
    if samples is None:
        samples = sorted(table.keys())
    
    # first, output the header/sample line
    yield "\t".join(['SEQUENCE'] + samples)
    
    # loop over rows
    for seq_label, sequence in lab_seq: 
        yield "\t".join([sequence] + [str(util_index.counts(table, sample, seq_label)) for sample in samples])


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('index', help='input index file')
    parser.add_argument('fasta', help='fasta file')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    parser.add_argument('--samples', default=None, help='samples in order in the first field (e.g., a barcode file)')
    
    args = parser.parse_args()
    
    # read in the index file
    with open(args.index) as f:
        table = sparse_seq_count_table(f)
    
    # read in the fasta file
    with open(args.fasta) as f:
        lab_seq = trimmed_fasta_entries(f)
    
    # read in the sample order, if given
    if args.samples is None:
        samples = None
    else:
        with open(args.samples) as f:
            samples = parse_sample_lines(f)
            
    for line in seq_table(table, lab_seq, samples):
        args.output.write("%s\n" % line)
