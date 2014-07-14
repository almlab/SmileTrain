#!/usr/bin/env python

'''
Create a sequence table from the index file. Option to use the dereplicated fasta
for ordering by abundance.
'''

import re, sys, argparse
import util, util_index


# adapted from uc2otus
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

# new
def full_seq_id_to_short_id(full_id):
    '''
    Convert usearch-style fasta id to just the sequence id.
    E.g., "seq123;size=4" to "seq123"
    '''
    
    return re.match("(seq\d+);", full_id).group(1)

def fasta_to_short_seq_ids(fasta_lines):
    '''Pull out the short (seq123) ids from a fasta file with long (seq123;size=400) ids'''
    return [full_seq_id_to_short_id(i) for i, s in util.fasta_entries(fasta_lines)]
    
# adapted
def seq_table(table, seqs=None, samples=None):
    '''
    Output lines of a seq table.
    
    table : dictionary of dictionaries
        {sample => {seq => abundance, ...}, ...}
    ordered_fasta : iterator of strings (default None)
        lines from the dereplicated fasta file
    samples : list or iterator of strings (default None)
        sample names in column order; or just make a new sorted list
        
    yields : strings
        lines in the otu table
    '''
    
    # make our own samples list if necessary
    if samples is None:
        samples = sorted(table.keys())
        
    # and our own seqs list
    if seqs is None:
        # concatenate the lists of keys
        seqs = []
        for sample in table:
            seqs += table[sample].keys()
            
        # remove duplicates and sort
        seqs = sorted(list(set(seqs)))
    
    # first, output the header/sample line
    yield "\t".join(['SEQ_ID'] + samples)
    
    # loop over rows
    for seq in seqs: 
        yield "\t".join([seq] + [str(util_index.counts(table, sample, seq)) for sample in samples])


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('index', help='input index file')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    parser.add_argument('--samples', default=None, help='samples in order in the first field (e.g., a barcode file)')
    parser.add_argument('--fasta', default=None, help='fasta file with sequences in order')
    args = parser.parse_args()
        
    with open(args.index) as f:
        table = sparse_seq_count_table(f)
        
    if args.samples is None:
        samples = None
    else:
        with open(args.samples) as f:
            samples = parse_sample_lines(f)
    
    if args.fasta is None:
        seqs = None
    else:
        with open(args.fasta) as f:
            seqs = fasta_to_short_seq_ids(f)
            
    for line in seq_table(table, seqs, samples):
        args.output.write("%s\n" % line)
