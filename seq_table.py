#!/usr/bin/env python

'''
Create a sequence file using the original and dereplicated fastas.

Load all the dereplicated sequences and sequence IDs into memory. Then read through the
original fasta. If a sequence is in the dereplicated fasta, then keep track of the
abundance of the sequence in that sample.
'''

import sys, argparse, re
import util, util_index
    
def fasta_to_table_and_abund(lines):
    '''
    Count the abundance of each sequence in each sample.
    
    lines : list or iterator of strings
        lines in the big fasta file
        
    returns : dict of dicts
        {sequence => {samples => counts}, ...}
    '''
    
    table = {}
    abund = {}
    for sid, seq in util.fasta_entries(lines):
        sample = util_index.sid_to_sample(sid)
        
        if seq in abund:
            abund[seq] += 1
        else:
            abund[seq] = 1
        
        if seq in table:
            if sample in table[seq]:
                table[seq][sample] += 1
            else:
                table[seq][sample] = 1
        else:
            table[seq] = {sample: 1}
            
    return table, abund

def table_to_samples(table):
    '''get sorted list of sample names'''
    samples = []
    for seq in table:
        for sample in table[seq]:
            if sample not in samples:
                samples.append(sample)
    
    samples = sorted(samples)
    return samples

def table_lines(fasta, min_counts, samples=None):
    '''fasta lines to seq table lines'''
    table, abund = fasta_to_table_and_abund(fasta)
    
    if samples is None:
        samples = table_to_samples(table)
    
    # get the sequences in abundance order
    all_seqs = sorted(table.keys(), key=lambda seq: abund[seq], reverse=True)
    
    # throw out sequences below minimum
    seqs = [seq for seq in all_seqs if abund[seq] >= min_counts]
    
    # write the header
    yield "sequence\t" + "\t".join(samples) + "\t" + "abundance"
    
    for seq in seqs:
        counts = [table[seq].get(sample, 0) for sample in samples]
        yield seq + "\t" + "\t".join([str(x) for x in counts]) + "\t" + str(abund[seq])


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('-s', '--samples', default=None, help='samples list (default: sorted names from fasta)')
    parser.add_argument('-m', '--minimum_counts', type=int, default=2, help='minimum times a sequence is included, otherwise it gets thrown out (default: 2)')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    args = parser.parse_args()
    
    if args.samples is not None:
        with open(args.samples) as f:
            samples = [line.strip() for line in f]
    else:
        samples = None
        
    with open(args.fasta) as f:
        for line in table_lines(f, args.minimum_counts, samples):
            print line