#!/usr/bin/env python

'''
Create a sequence file using the original and dereplicated fastas.

Load all the dereplicated sequences and sequence IDs into memory. Then read through the
original fasta. If a sequence is in the dereplicated fasta, then keep track of the
abundance of the sequence in that sample.
'''

import sys, argparse, re
from Bio import SeqIO
import util, util_index

def fasta_to_dict(fasta):
    '''
    Create a hash {sequence => ID}

    fasta : fasta fh or fn
        lines in the fasta

    returns : dict
        {sequence => ID}
    '''

    names = {str(record.seq): record.id for record in SeqIO.parse(fasta, 'fasta')}
    return names

def fasta_to_table_and_abund(fasta):
    '''
    Count the abundance of each sequence in each sample.
    
    fasta : fasta fh or fn
        lines in the big fasta file
        
    returns : dict of dicts
        {sequence => {samples => counts}, ...}
    '''
    
    table = {}
    abund = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        sample = util_index.sid_to_sample(record.id)
        seq = str(record.seq)
        
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

def table_lines(fasta, derep, min_counts, samples=None):
    '''fasta lines to seq table lines'''

    # read in the dereplicated fasta
    names = fasta_to_dict(derep)

    table, abund = fasta_to_table_and_abund(fasta)

    # check that all the sequences in the big file are in the dereplicated and vice versa
    for seq in table:
        if seq not in names:
            raise RuntimeError("sequence found in full but no dereplicated fasta:\n%s" %(seq))

    for seq_id, seq in names.items():
        if seq not in table:
            raise RuntimeError("sequence %s found in dereplicated but not full fasta" %(seq))
    
    if samples is None:
        samples = table_to_samples(table)
    
    # get the sequences in abundance order
    all_seqs = sorted(table.keys(), key=lambda seq: abund[seq], reverse=True)
    
    # throw out sequences below minimum
    seqs = [seq for seq in all_seqs if abund[seq] >= min_counts]
    
    # write the header
    yield "sequence_id\t" + "\t".join(samples)
    
    for seq in seqs:
        counts = [table[seq].get(sample, 0) for sample in samples]
        yield names[seq] + "\t" + "\t".join([str(x) for x in counts])


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='full fasta file')
    parser.add_argument('derep', help='dereplicated fasta file')
    parser.add_argument('-s', '--samples', default=None, help='samples list (samples in first field; default: sorted names from fasta)')
    parser.add_argument('-m', '--minimum_counts', type=int, default=2, help='minimum times a sequence is included, otherwise it gets thrown out (default: 2)')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    args = parser.parse_args()
    
    if args.samples is not None:
        with open(args.samples) as f:
            samples = [line.split()[0] for line in f]
    else:
        samples = None
        
    for line in table_lines(args.fasta, args.derep, args.minimum_counts, samples):
        args.output.write(line + "\n")