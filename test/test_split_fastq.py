#!/usr/bin/env python

'''
for split_fastq.py
'''

from SmileTrain.test import fake_fh
import unittest, StringIO
from Bio import SeqIO, Seq

from SmileTrain import split_fastq


class TestOutputFilenames(unittest.TestCase):
    def test_correct(self):
        self.assertEqual(split_fastq.output_filenames('foo.fasta', 3), ['foo.fasta.0', 'foo.fasta.1', 'foo.fasta.2'])


class TestSplitFastqEntries(unittest.TestCase):
    def test_correct(self):
        in_fh = fake_fh(['@foo', 'AAA', '+foo', '###', '@bar', 'CCC', '+bar', '"""', '@baz', 'TTT', '+baz', '$$$', '@poo', 'GGG', '+poo', '==='])
        outs = [fake_fh() for x in range(3)]
        
        split_fastq.split_fastq_entries(in_fh, outs)
        conts = [out.getvalue() for out in outs]
        
        exp1 = "\n".join(['@foo', 'AAA', '+', '###', '@poo', 'GGG', '+', '===']) + "\n"
        exp2 = "\n".join(['@bar', 'CCC', '+', '"""']) + "\n"
        exp3 = "\n".join(['@baz', 'TTT', '+', '$$$']) + "\n"
        
        self.assertEqual(conts, [exp1, exp2, exp3])


if __name__ == '__main__':
    unittest.main(verbosity=2)
