#!/usr/bin/env python

'''
unit tests for split_fasta.py
'''

from SmileTrain.test import fake_fh
import unittest
from SmileTrain import split_fasta

class TestSplitFastaEntries(unittest.TestCase):
    def setUp(self):
        self.fh = fake_fh(">foo\nAAA\n>bar\nCCC\n>baz\nTTT\n>poo\nGGG\n")
    
    def test_correct(self):
        outs = [fake_fh() for x in range(3)]
        split_fasta.split_fasta_entries(self.fh, outs)
        conts = [out.getvalue() for out in outs]
        self.assertEqual(conts, [">foo\nAAA\n>poo\nGGG\n", ">bar\nCCC\n", ">baz\nTTT\n"])
    
    def test_correct_by_hash(self):
        outs = [fake_fh() for x in range(2)]
        split_fasta.split_fasta_entries(self.fh, outs, by_hash=True)
        conts = [out.getvalue() for out in outs]
        self.assertEqual(conts, [">foo\nAAA\n>bar\nCCC\n>poo\nGGG\n", ">baz\nTTT\n"])


if __name__ == '__main__':
    unittest.main(verbosity=2)
