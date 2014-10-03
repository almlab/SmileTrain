import unittest
from SmileTrain import util_index

class TestParseIndexLine(unittest.TestCase):
    '''tests for function parse_index_line'''
    
    def test_correct(self):
        '''should split and recast index line'''
        line = 'donor1  seq5    100'
        self.assertEqual(util_index.parse_index_line(line), ['donor1', 'seq5', 100])
     
   
class TestCounts(unittest.TestCase):
    '''tests for function counts'''

    def setUp(self):
        self.table = {'sample1': {'otu1': 1, 'otu2': 2}, 'sample2': {'otu1': 2, 'otu2': 5}}
    
    def test_correct_present(self):
        '''should parse the counts table structure'''
        self.assertEqual(util_index.counts(self.table, 'sample2', 'otu2'), 5)

    def test_correct_absent(self):
        '''should give 0 if can't find sample and otu'''
        self.assertEqual(util_index.counts(self.table, 'sample2', 'otu3'), 0)
    

class TestParseSeqSid(unittest.TestCase):
    def test_correct(self):
        self.assertEqual(util_index.parse_seq_sid('seq44;counts=12'), 'seq44')