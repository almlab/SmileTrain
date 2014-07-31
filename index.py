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
from Bio import SeqIO
import util, util_index


def parse_derep_fasta(fasta):
    '''create a hash {sequence => ID} from fasta filename or filehandle'''
    return {str(record.seq): util_index.parse_seq_id(record.id) for record in SeqIO.parse(fasta, 'fasta')}

def sid_to_sample(sid):
    '''sample=donor1;400 -> donor1'''
    
    m = re.match('sample=(.+);\d+', sid)
    if m is None:
        raise RuntimeError("fasta at line did not parse: %s" % sid)
    else:
        return m.group(1)
    
def parse_full_fasta(fasta, seq_sid):
    '''
    Count the abundance of each sequence in each sample. Ignore sequences that are not
    in the known ID list.
    
    fasta : filehandle or filename
        input fasta
    seq_sid : dictionary
        {sequence => sequence ID}
        
    returns : dict
        {(sample, sequence ID) => abundance}
    '''
    
    abund = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        sample = sid_to_sample(record.id)
        
        if record.seq in seq_sid:
            seq_id = seq_sid[str(record.seq)]
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
    
    seq_sid = parse_derep_fasta(args.derep)
    abundances = parse_full_fasta(args.orig, seq_sid)
                
    for line in index_lines(abundances):
        args.output.write(line + "\n")