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

If if by_sample=True, then keep a list of where the counts came from to produce an
index file with lines
    sample_name sequence_id (# of times that seq appears in that sample)
'''

import sys, argparse, re
import util

class Dereplicator():
    def __init__(self, fasta_entries, minimum_counts, by_sample):
        '''
        fasta_entries : list or iterator
            [sid (no >), sequence] pairs
        minimum_counts : int
            minimum number of counts to be included in output
        by_sample : bool
            dereplicate only within samples?
        '''
        
        self.input_fasta_entries = fasta_entries
        self.minimum_counts = minimum_counts
        self.by_sample = by_sample
        
        # process the data
        self.dereplicate()
        
        # sort out the abundant sequences
        self.sort_abundant_sequences()

    def sid_to_sample(self, sid):
        '''sample=donor1;400 -> donor1'''
        
        m = re.match('sample=(.+);\d+', sid)
        if m is None:
            raise RuntimeError("fasta at line did not parse: %s" % sid)
        else:
            return m.group(1)
        
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
            
        If dereplicating by sample, the abundance of each sequence in each smaple is kept
            sample_counts = {('donor1', 'seq1' => 40}
        '''
        
        self.abundances = {}
        self.seq_ids = {}
        new_seq_ids = self.iter_seq_ids()
        
        if self.by_sample:
            self.sample_counts = {}
        
        # keep track of the highest sequence index used
        for sid, seq in self.input_fasta_entries:    
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
                
            # if dereplicating by sample, keep track of the abundance of this sequence in
            # its own sample
            if self.by_sample:
                sample = self.sid_to_sample(sid)
                key = (sample, seq_id)
                if key not in self.sample_counts:
                    self.sample_counts[key] = 1
                else:
                    self.sample_counts[key] += 1
                    
    def sort_abundant_sequences(self):
        '''get the abundant sequences in order of their abundance'''
        
        sorted_seqs = sorted(self.seq_ids.keys(), key=lambda seq: self.abundances[seq], reverse=True)
        self.filtered_abundant_sequences = [seq for seq in sorted_seqs if self.abundances[seq] >= self.minimum_counts]
        self.filtered_abundant_ids = [self.seq_ids[seq] for seq in self.filtered_abundant_sequences]
    
    def seq_to_entry(self, seq):
        '''seq_id, abundance -> >seq_id;counts=abundance\nACGT'''
        
        return ">%s;counts=%d\n%s" %(self.seq_ids[seq], self.abundances[seq], seq)
    
    def new_fasta_entries(self):
        '''yield the list of fasta lines in abundance order'''
        
        for seq in self.filtered_abundant_sequences:
            yield self.seq_to_entry(seq)
            
    def sample_index_entries(self):
        '''yield the list of lines in the index file'''
        
        for (sample, seq_id), abundance in self.sample_counts.items():
            if seq_id in self.filtered_abundant_ids:
                yield "%s\t%s\t%d" %(sample, seq_id, abundance)

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input fasta file')
    parser.add_argument('output', help='output fasta file')
    parser.add_argument('--index', type=str, help='if specified, create a sample index file')
    parser.add_argument('-m', '--minimum_counts', type=int, default=2, help='minimum times a sequence is included, otherwise it gets thrown out (default: 2)')
    args = parser.parse_args()
    
    if args.index:
        make_index = True
    else:
        make_index = False

    with open(args.input, 'r') as f:
        derep = Dereplicator(util.fasta_entries(f), args.minimum_counts, make_index)
        
    with open(args.output, 'w') as o:
        for entry in derep.new_fasta_entries():
            o.write("%s\n" %(entry))
            
    if make_index:
        with open(args.index, 'w') as o:
            for entry in derep.sample_index_entries():
                o.write("%s\n" %(entry))