import unittest, tempfile, subprocess, os, shutil
from SmileTrain import uc2denovo, util

class TestMissedLabels(unittest.TestCase):
    '''tests for missing_labels'''
    def test_correct(self):
        '''should parse uc file and get missing labels'''
        uc_lines = ["H       43107   253     100.0   +       0       0       526I253M617I    seq1;counts=617550      181719", "H       42869   252     100.0   +       0       0       526I252M625I    seq24;counts=296908     182945", "N       *       *       *       .       *       *       *       seq12996;counts=2       *"]
        labels = uc2denovo.missed_labels(uc_lines)
        self.assertEqual(labels, ['seq12996;counts=2'])
      
        
class TestMatchingFastaEntries(unittest.TestCase):
    '''tests for matching_fasta_entries'''
    def test_correct(self):
        '''should collect only the matching entries'''
        old_fasta_lines = [">a", "AAA", ">b", "TTT", ">c", "GGG"]
        labels = ["a", "c"]
        expected_new_fasta = ">a\nAAA\n>c\nGGG\n"
        
        fastas = util.fasta_entries(old_fasta_lines)
        new_fasta = uc2denovo.matching_fasta_entries(labels, fastas)
        
        self.assertEqual(new_fasta, expected_new_fasta)
        