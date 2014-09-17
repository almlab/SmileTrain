#!/usr/bin/env python

'''
Split an input fasta file into k files.

usage:
  python split_fasta.py in.fst 3

creates 3 files:
  in.fst.1
  in.fst.2
  in.fst.3
  
'''

import itertools, os.path, sys, argparse
from Bio import SeqIO
import util

def output_filenames(input_filename, k):
    '''destination filenames foo.fastq.0, etc.'''
    return ['%s.%d' % (input_filename, i) for i in range(k)]

def split_fasta_entries(fasta, fhs, by_hash=False):
    '''
    Send entries in the input to filenames, cycling over each filename.
    
    fasta : filehandle or filename
        input fasta file
    fhs : list or iterator of filehandles 
        output filehandles
    by_hash : bool (default false)
        split based on hash value rather than just cycling
        
    returns : nothing
    '''
    
    if by_hash:
        # make the function for determining which bin the sequences fall in
        h = lambda seq: hash(seq) % len(fhs)
        
        # pick the filehandle bashed on the sequence's hash, then write
        for record in SeqIO.parse(fasta, 'fasta'):
            SeqIO.write(record, fhs[h(str(record.seq))], 'fasta')
    else:
        fh_cycler = itertools.cycle(fhs)
        
        for record, fh in itertools.izip(SeqIO.parse(fasta, 'fasta'), fh_cycler):
            SeqIO.write(record, fh, 'fasta')

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split a fasta file foo.fasta into multiple fastq files foo.fasta.0, foo.fasta.1, etc.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='input fasta')
    parser.add_argument('n_files', type=int, help='number of split files to output')
    parser.add_argument('-s', '--hash', action='store_true', help='split by hash')
    args = parser.parse_args()
    
    filenames = output_filenames(args.fasta, args.n_files)
    
    util.check_for_collisions(filenames)
    
    if len(filenames) == 1:
        # just copy the file
        shutil.copy(args.fasta, filenames[0])
    else:
        # split the file entry by entry
        split_fasta_entries(args.fasta, [open(f, 'w') for f in filenames], by_hash=args.hash)