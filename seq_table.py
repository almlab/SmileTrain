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

class SeqTableWriter:
    def __init__(self, fasta, derep, output, samples=None, min_counts=0, assert_same_seqs=False, run=True):
        self.fasta = fasta
        self.derep = derep
        self.output = output
        self.samples = samples
        self.min_counts = min_counts
        self.assert_same_seqs = assert_same_seqs

        if run:
            self.run()

    def run(self):
        self.derep_names = self.fasta_to_dict(self.derep)
        self.table, self.abund = self.fasta_to_abund(self.fasta, self.derep_names, self.assert_same_seqs)

        if self.samples is None:
            self.samples = self.table_to_samples(self.table)

        self.write_table(self.table, self.abund, self.samples, self.min_counts, self.output)

    @staticmethod
    def fasta_to_dict(fasta):
        '''
        Create a hash {sequence => ID}

        fasta : fasta fh or fn
            lines in the fasta

        returns : dict
            {sequence => ID}
        '''

        names = {str(record.seq): util_index.parse_seq_sid(record.id) for record in SeqIO.parse(fasta, 'fasta')}
        return names

    @staticmethod
    def fasta_to_abund(fasta, names, assert_same=False):
        '''
        Count the abundance of each sequence in each sample.
        
        fasta : fasta fh or fn
            lines in the big fasta file
        names : dict
            {seq => name}
        assert_same : bool
            if true, make sure each seq in the fasta is in names
            
        returns : dict of dicts
            {name => {samples => counts}, ...}
        '''
    
        table = {}
        abund = {}
        for record in SeqIO.parse(fasta, 'fasta'):
            sample = util_index.sid_to_sample(record.id)
            seq = str(record.seq)

            if seq in names:
                name = names[seq]
            else:
                if assert_same:
                    raise RuntimeError("sequence %s found in fasta but not dereplicated fasta" %(seq))

            if name in abund:
                abund[name] += 1
            else:
                abund[name] = 1
            
            if name in table:
                if sample in table[name]:
                    table[name][sample] += 1
                else:
                    table[name][sample] = 1
            else:
                table[name] = {sample: 1}
                
        return table, abund

    @staticmethod
    def table_to_samples(table):
        '''get sorted list of sample names'''
        samples = []
        for seq in table:
            for sample in table[seq]:
                if sample not in samples:
                    samples.append(sample)
        
        samples = sorted(samples)
        return samples

    @staticmethod
    def write_table(table, abund, samples, min_counts, output):
        '''fasta lines to seq table lines'''

        # get the seq ids in abundance order
        all_seq_ids = sorted(table.keys(), key=lambda i: abund[i], reverse=True)
        
        # throw out sequences below minimum
        seq_ids = [i for i in all_seq_ids if abund[i] >= min_counts]
        
        # write the header
        output.write("sequence_id\t" + "\t".join(samples) + "\n")
        
        for i in seq_ids:
            counts = [table[i].get(sample, 0) for sample in samples]
            output.write(i + "\t" + "\t".join([str(x) for x in counts]) + "\n")


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Make a sequence table (i.e., 100% identity OTUs)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='full fasta file')
    parser.add_argument('derep', help='dereplicated fasta file')
    parser.add_argument('-s', '--samples', default=None, help='samples list (samples in first field; default: sorted names from fasta)')
    parser.add_argument('-m', '--min_counts', type=int, default=0, help='minimum times a sequence is included, otherwise it gets thrown out')
    parser.add_argument('-a', '--assert_same_seqs', action='store_true', help='assert that every seq in fasta is in dereplicated fasta?')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='output file')
    args = parser.parse_args()

    opts = vars(args)
    SeqTableWriter(**opts)