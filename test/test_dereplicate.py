#!/usr/bin/env python

'''
unit tests for derep_fulllength.py
'''

from SmileTrain.test import fake_fh
import unittest
from SmileTrain import derep_fulllength

class TestDereplicate(unittest.TestCase):
    '''tests for dereplication without samples'''
    def setUp(self):
        fasta = fake_fh(['>foo', 'AA', '>bar', 'CC', '>baz', 'AA', '>blag', 'AA', '>flog', 'TT', '>blob', 'TT'])
        minimum_counts = 2
        self.derep = derep_fulllength.Dereplicator(fasta, minimum_counts)
        
    def test_new_seq_id(self):
        '''should give increasing sequence ids'''
        it = self.derep.iter_seq_ids()
        self.assertEqual(it.next(), "seq0")
        self.assertEqual(it.next(), "seq1")
        self.assertEqual(it.next(), "seq2")
        
    def test_dereplicate_abundances(self):
        '''should properly count abundances'''
        self.assertEqual(self.derep.abundances, {'AA': 3, 'CC': 1, 'TT': 2})
        
    def test_dereplicate_seq_ids(self):
        '''should properly organize names'''
        self.assertEqual(self.derep.seq_ids, {'AA': 'seq0', 'CC': 'seq1', 'TT': 'seq2'})
        
    def test_sort_abundant_sequences(self):
        '''should properly short abundances'''
        self.assertEqual(self.derep.filtered_abundant_sequences, ['AA', 'TT'])
    
    def test_fasta_entries(self):
        '''should give abundance-sorted entries'''
        fe = self.derep.new_fasta_entries()
        record1 = fe.next()
        record2 = fe.next()
        self.assertEqual([record1.id, str(record1.seq), record2.id, str(record2.seq)], ['seq0;counts=3', 'AA', 'seq2;counts=2', 'TT'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
