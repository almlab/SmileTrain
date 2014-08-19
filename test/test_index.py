#!/usr/bin/env python

'''
unit tests for ???
'''

from SmileTrain.test import fake_fh
import unittest
from SmileTrain import index

class TestIndex(unittest.TestCase):
    '''tests for index-writing script'''
        
    def test_parse_derep_fasta(self):
        '''should make a dictionary of fasta lines'''
        fasta = fake_fh(['>seq0;counts=10', 'AAA', '>seq4;counts=23', 'TTT'])
        self.assertEqual(index.parse_derep_fasta(fasta), {'AAA': 'seq0', 'TTT': 'seq4'})
        
    def test_sid_to_sample(self):
        '''should extract sample from fasta line'''
        self.assertEqual(index.sid_to_sample('sample=donor1;444'), 'donor1')
        
    def test_parse_full_fasta(self):
        seq_sid = {'AAA': 'seq0', 'TTT': 'seq4'}
        fasta = fake_fh(['>sample=donor1;1', 'AAA', '>sample=donor1;2', 'AAA', '>sample=donor1;3', 'TTT', '>sample=donorT;1', 'TTT'])
        abund = index.parse_full_fasta(fasta, seq_sid)
        self.assertEqual(abund, {('donor1', 'seq0'): 2, ('donor1', 'seq4'): 1, ('donorT', 'seq4'): 1})



if __name__ == '__main__':
    unittest.main(verbosity=2)
