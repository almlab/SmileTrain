#!/usr/bin/env python

'''
Dereplicate a fasta file. Read through all the entries, finding reads with identical
sequences. Write out a new fasta file with IDs that keep track of the abundance of
each sequence.

input like
    >foo
    CATCAT
    >bar
    TATA
    >baz
    AGGAGG
    >lorem
    TATA
    >ipsum
    CATCAT
    >qualcum
    CATCAT

with a minimum number of reads
    2

becomes output like
    >seq1;size=3
    CATCAT
    >seq2;size=2
    TATA

The output file is sorted in order of decreasing abundance. Sequences with a minimum
number of counts are dropped.
'''

import sys, argparse, re, sys
from Bio import SeqIO, SeqRecord, Seq
import util

class Dereplicator():
    def __init__(self, fasta, minimum_counts):
        '''
        fasta : fasta filename or handle
            input fasta
        minimum_counts : int
            minimum number of counts to be included in output
        '''
        
        self.fasta = fasta
        self.minimum_counts = minimum_counts
        
        # process the data
        self.dereplicate()
        
        # sort out the abundant sequences
        self.sort_abundant_sequences()
        
    def iter_seq_ids(self):
        self.max_seq_i = 0
        while True:
            yield "seq%d" % self.max_seq_i
            self.max_seq_i += 1
        
    def dereplicate(self):
        '''
        Process the input fasta entries.
        
        Each unique sequence is given a unique id
            seq_ids = {'ACGT' => 'seq1'}
            
        Each sequence also has an abundance
            abundances = {'ACGT' => 100}
        '''
        
        self.abundances = {}
        self.seq_ids = {}
        new_seq_ids = self.iter_seq_ids()
        
        # keep track of the highest sequence index used
        for record in SeqIO.parse(self.fasta, 'fasta'):
            seq = str(record.seq)
            
            # if we haven't seen this sequence before, give it a new index
            if seq not in self.seq_ids:
                # this should be a new id
                new_id = new_seq_ids.next()
                assert(new_id not in self.seq_ids)
                
                self.seq_ids[seq] = new_id
                
            # get the id for this sequence
            seq_id = self.seq_ids[seq]
            
            # track the abundance of this sequence
            if seq not in self.abundances:
                self.abundances[seq] = 1
            else:
                self.abundances[seq] += 1
                    
    def sort_abundant_sequences(self):
        '''get the abundant sequences in order of their abundance'''
        
        sorted_seqs = sorted(self.seq_ids.keys(), key=lambda seq: self.abundances[seq], reverse=True)
        self.filtered_abundant_sequences = [seq for seq in sorted_seqs if self.abundances[seq] >= self.minimum_counts]
        self.filtered_abundant_ids = [self.seq_ids[seq] for seq in self.filtered_abundant_sequences]
    
    def seq_to_record(self, seq):
        '''seq_id, abundance -> >seq_id;counts=abundance\nACGT'''
        return SeqRecord.SeqRecord(Seq.Seq(seq), id="%s;counts=%d" %(self.seq_ids[seq], self.abundances[seq]), description='')
    
    def new_fasta_entries(self):
        '''yield the list of fasta lines in abundance order'''
        
        for seq in self.filtered_abundant_sequences:
            yield self.seq_to_record(seq)


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Make a new fasta, keeping only one entry for every unique sequence', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='input fasta')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='dereplicated fasta')
    parser.add_argument('-m', '--minimum_counts', type=int, default=2, help='minimum times a sequence is included, otherwise it gets thrown out')
    args = parser.parse_args()

    derep = Dereplicator(args.fasta, args.minimum_counts)
    
    SeqIO.write(derep.new_fasta_entries(), args.output, 'fasta')