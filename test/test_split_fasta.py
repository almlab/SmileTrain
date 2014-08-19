#!/usr/bin/env python

'''
unit tests for split_fasta.py
'''

from SmileTrain.test import fake_fh
import unittest
from SmileTrain import split_fasta

class TestSplitFastaEntries(unittest.TestCase):
    def setUp(self):
        self.fh = fake_fh(">foo\nAAA\n>bar\nCCC\n>baz\nTTT\npoo\nGGG\n")
    
    def test_correct(self):
        outs = [fake_fh() for x in range(3)]
        split_fasta.split_fasta_entries(self.fh, outs, 3)
    
    def test_correct_by_hash(self):
        raise RuntimeError("test not implemented")


if __name__ == '__main__':
    unittest.main(verbosity=2)
