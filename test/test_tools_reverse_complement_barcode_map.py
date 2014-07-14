import unittest, tempfile, subprocess, os, shutil
import Bio.Seq
from SmileTrain.tools import reverse_complement_barcode_map as rcbm

class TestReverseComplementMap(unittest.TestCase):
    '''tests for function reverse_complement_map'''
    
    def test_correct(self):
        '''should keep the sample name and reverse the barcode'''
        old_map = [['sample1', 'AATT'], ['sample2', 'CCTT']]
        new_map = [['sample1', 'AATT'], ['sample2', 'AAGG']]
        self.assertEqual(new_map, rcbm.reverse_complement_map(old_map))


if __name__ == '__main__':
    unittest.main(verbosity=2)
