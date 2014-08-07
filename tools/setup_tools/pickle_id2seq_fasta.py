#!/usr/bin/env python

'''
Pickle fasta file into a dictionary {id => sequence}.
'''

import argparse, cPickle as pickle
from Bio import SeqIO


def fasta_file_to_dict(fh):
    '''fasta file to {id => sequence}'''
    return {record.id: str(record.seq) for record in SeqIO.parse(fh, 'fasta')}


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta')
    parser.add_argument('-o', '--output', default='fasta.pkl', help='output pickle file (default: fasta.pkl')
    args = parser.parse_args()
    
    d = fasta_file_to_dict(input.fasta)
    
    # write the dictionary as a binary pickle
    with open(args.output, 'wb') as f:
        pickle.dump(d, f)