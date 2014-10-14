import unittest, tempfile, subprocess, os, shutil
from SmileTrain import util, remove_primers, derep_fulllength, check_fastq_format, convert_fastq, map_barcodes, derep_fulllength, uc2otus, index
from SmileTrain.test import fake_fh


class TestBarcodeDictionary(unittest.TestCase):
    '''tests for the barcode dictionary'''
    def test_correct(self):
        '''should parse a tab-separated mapping file'''
        mapping_lines = ['donor1\tACGT', 'donor2\tTACA']
        d = map_barcodes.barcode_file_to_dictionary(mapping_lines)
        self.assertEqual(d, {'ACGT': 'donor1', 'TACA': 'donor2'})

    def test_incomplete_line(self):
        '''should raise an error for an incomplete line'''
        mapping_lines = ['donor1\tACGT', 'donor2\t']
        self.assertRaises(RuntimeError, map_barcodes.barcode_file_to_dictionary, mapping_lines)


class TestBestBarcode(unittest.TestCase):
    def test_correct(self):
        '''should find the best barcode out of a list'''
        barcodes = ['AAAATT', 'TACACC', 'ACGTAA']
        target = 'TACTCC'
        n_mismatches, best = map_barcodes.best_barcode_match(barcodes, target)
        self.assertEqual(n_mismatches, 1)
        self.assertEqual(best, 'TACACC')


class TestRename(unittest.TestCase):
    def test_correct(self):
        '''should properly rename samples'''
        d = {'ACGT': 'donor1', 'TACA': 'donor2'}
        fastq = fake_fh(['@foo#ACGT/1', 'AAA', '+', 'AAA', '@bar#TACA/1', 'CCC', '+', 'BBB'])
        it = map_barcodes.renamed_fastq_records(fastq, d, 1)
        
        record = it.next()
        self.assertEqual(str(record.seq), 'AAA')
        self.assertEqual(record.id, 'sample=donor1;1/1')
        record = it.next()
        self.assertEqual(str(record.seq), 'CCC')
        self.assertEqual(record.id, 'sample=donor2;1/1')
        
    def test_no_barcode(self):
        '''should raise error when no barcode in the @ line'''
        d = {'ACGT': 'donor1', 'TACA': 'donor2'}
        fastq = fake_fh(['@foo#ACGT/1', 'AAA', '+', 'AAA', '@bar_bad_format/1', 'CCC', '+', 'BBB'])
        it = map_barcodes.renamed_fastq_records(fastq, d, 1)
        
        # the first entry is OK, but the second should complain
        it.next()
        self.assertRaises(RuntimeError, it.next)


if __name__ == '__main__':
    unittest.main(verbosity=2)
