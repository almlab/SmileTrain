#!/usr/bin/env python

'''
unit tests for check_fastq_format.py
'''

from SmileTrain.test import fake_fh
import unittest
from SmileTrain import check_fastq_format
from Bio import SeqIO

tmp_dir = 'tmp'

class FastqTest(unittest.TestCase):
    '''tests for functions and scripts that do fastq manipulations'''
    def setUp(self):
        self.good_fastq_content = "@foo\nAAA\n+foo\n!!!\n@bar\nCCC\n+bar\n###"
        self.fastq13 = """@lolapolooza:1234#ACGT/1\nAATTAAGTCAAATTTGGCCTGGCCCAGTGTCCAATGTTGT\n+\nABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh\n"""
        self.fastq18 = """@lolapolooza:1234#ACGT/1\nAATTAAGTCAAATTTGGCCTGGCCCAGTGTCCAATGTTGT\n+\n"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI\n"""
        
        self.fastq13_fh = fake_fh(self.fastq13 + "\n")
        self.fastq13_record = SeqIO.read(self.fastq13_fh, 'fastq')
        
        self.fastq18_fh = fake_fh(self.fastq18 + "\n")
        self.fastq18_record = SeqIO.read(self.fastq18_fh, 'fastq')
        
        fh = fake_fh('''@lol:1234#ACGT/1\nAAAA\n+\nAh"J\n''')
        self.fastq_ambig_record = SeqIO.read(fh, 'fastq')


class TestParseFastqRecordID(FastqTest):
    def test_correct(self):
        '''fastq_at_line_to_id should trim the starting @ and trailing /1 or /2'''
        self.assertEqual(check_fastq_format.parse_fastq_record_id(self.fastq13_record), "lolapolooza:1234#ACGT")


class TestFastqRecordFormat(FastqTest):
    def test_illumina13_id(self):
        '''should find illumina13 format'''
        self.assertEqual(check_fastq_format.fastq_record_format(self.fastq13_record), 'illumina13')
        
    def test_illumina18_id(self):
        '''should find illumina18 format'''
        self.assertEqual(check_fastq_format.fastq_record_format(self.fastq18_record), 'illumina18')
        
    def test_ambiguous_illumina_id(self):
        '''should find ambiguous illumina format'''
        self.assertEqual(check_fastq_format.fastq_record_format(self.fastq_ambig_record), 'ambiguous')
        
        
class TestCheckIlluminaFormat(FastqTest):
    def test_correct(self):
        '''should identify written file as correct format'''
        with open('%s/good.fastq' % tmp_dir, 'w') as f:
            f.write(self.fastq13)
        
        check_fastq_format.check_illumina_format(['%s/good.fastq' % tmp_dir], 'illumina13')
        
    def test_incorrect(self):
        '''should identify written file with incorrect format'''
        with open('%s/bad.fastq' % tmp_dir, 'w') as f:
            f.write(self.fastq18)
        
        self.assertRaises(RuntimeError, check_fastq_format.check_illumina_format, ['%s/bad.fastq' % tmp_dir], 'illumina13')
            
    def test_format_conversion(self):
        '''should convert Illumina 1.4-1.7 to our mixed format'''
        
        record18 = convert_fastq.convert_record_illumina13_to_18(self.fastq13_record)
        fh = fake_fh()
        SeqIO.write(record18, fh, 'fastq')
        entry18 = fh.getvalue()
            
        self.assertEqual(entry18, self.fastq18)


if __name__ == '__main__':
    unittest.main(verbosity=2)
