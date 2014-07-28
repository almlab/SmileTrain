#!/usr/bin/env python

'''
Create an index file using the original and dereplicated fastas.

Load all the dereplicated sequences and sequence IDs into memory. Then read through the
original fasta. If a sequence is in the dereplicated fasta, then keep track of the
abundance of the sequence in that sample.

Produce an index file with lines
    sample_name sequence_id (# of times that seq appears in that sample)
'''

import sys, argparse, re
import util, util_index


def parse_derep_fasta(lines):
    '''create a hash {sequence => ID}'''
    return {seq: util_index.parse_seq_sid(sid) for sid, seq in util.fasta_entries(lines)}

def sid_to_sample(sid):
    '''sample=donor1;400 -> donor1'''
    
    m = re.match('sample=(.+);\d+', sid)
    if m is None:
        raise RuntimeError("fasta at line did not parse: %s" % sid)
    else:
        return m.group(1)
    
def parse_full_fasta(lines, seq_sid):
    '''
    Count the abundance of each sequence in each sample. Ignore sequences that are not
    in the known ID list.
    
    lines : list or iterator of strings
        lines in the big fasta file
    seq_sid : dictionary
        {sequence => sequence ID}
        
    returns : dict
        {(sample, sequence ID) => abundance}
    '''
    
    abund = {}
    for sid, seq in util.fasta_entries(lines):
        sample = sid_to_sample(sid)
        
        if seq in seq_sid:
            seq_id = seq_sid[seq]
            key = (sample, seq_id)
            
            if key in abund:
                abund[key] += 1
            else:
                abund[key] = 1
                
    return abund
                
def index_lines(abundances):
    '''{(sample, ID) => abundance, ...} -> tab-separated sample, ID, abundance'''
    
    for (sample, sid), abundance in abundances.items():
        yield "\t".join([sample, sid, str(abundance)])


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('orig', help='original fasta file')
    parser.add_argument('derep', help='dereplicated fasta file')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    args = parser.parse_args()
    
    with open(args.derep) as f:
        seq_sid = parse_derep_fasta(f)
                
    with open(args.orig) as f:
        abundances = parse_full_fasta(f, seq_sid)
                
    for line in index_lines(abundances):
        args.output.write(line + "\n")
