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
    
    def test_correct(self):
        '''should parse the counts table structure'''
        raise RuntimeError("test not implemented")