#!/usr/bin/env python

'''
unit tests for convert_fastq.py
'''

from SmileTrain.test import fake_fh
import unittest
from SmileTrain import convert_fastq
from Bio import SeqIO


class FastqTest(unittest.TestCase):
    '''tests for functions and scripts that do fastq manipulations'''
    def setUp(self):
        self.fastq13 = """@lolapolooza:1234#ACGT/1\nAATTAAGTCAAATTTGGCCTGGCCCAGTGTCCAATGTTGT\n+\nABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh\n"""
        self.fastq18 = """@lolapolooza:1234#ACGT/1\nAATTAAGTCAAATTTGGCCTGGCCCAGTGTCCAATGTTGT\n+\n"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI\n"""
        
        self.fastq13_fh = fake_fh(self.fastq13 + "\n")
        self.fastq13_record = SeqIO.read(self.fastq13_fh, 'fastq')


class TestConvertRecord(FastqTest):
    def test_correct(self):
        '''should convert Illumina 1.4-1.7 to our mixed format'''
        
        record18 = convert_fastq.convert_record_illumina13_to_18(self.fastq13_record)
        fh = fake_fh()
        SeqIO.write(record18, fh, 'fastq')
        entry18 = fh.getvalue()
            
        self.assertEqual(entry18, self.fastq18)


if __name__ == '__main__':
    unittest.main(verbosity=2)
