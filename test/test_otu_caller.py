#!/usr/bin/env python

'''
Unit tests for otu_caller.py
'''

from SmileTrain.test import fake_fh
import unittest, tempfile, subprocess, os, shutil, StringIO
from Bio import SeqIO, Seq

from SmileTrain import util, remove_primers, derep_fulllength, check_fastq_format, convert_fastq, map_barcodes, derep_fulllength, uc2otus, index, otu_caller

tmp_dir = 'tmp'

class TestWithFiles(unittest.TestCase):
    '''tests that need to read and write files'''
    
    def setUp(self):
        os.mkdir(tmp_dir)
    
    def tearDown(self):
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)

            
class TestFileChecks(TestWithFiles):
    '''tests for utilities that make queries about files'''
    
    def setUp(self):
        os.mkdir(tmp_dir)
        empty_fh, self.empty_fn = tempfile.mkstemp(dir=tmp_dir)
        
        full_fh, self.full_fn = tempfile.mkstemp(dir=tmp_dir)
        os.write(full_fh, "hello world")
        os.close(full_fh)
        
        no_fh, self.no_fn = tempfile.mkstemp(dir=tmp_dir)
        os.close(no_fh)
        os.unlink(self.no_fn)

        self.sub = otu_caller.Submitter(method='local')
        
    def test_check_for_existence_empty(self):
        '''should identify empty files as existing'''
        self.sub.check_for_existence(self.empty_fn)
    
    def test_check_for_existence_full(self):
        '''should identify full file as existing'''
        self.sub.check_for_existence(self.full_fn)
        
    def test_check_for_existence_no(self):
        '''should identify deleted file as missing'''
        self.assertRaises(RuntimeError, self.sub.check_for_existence, self.no_fn)
            
    def test_check_for_collision_empty(self):
        '''should identify empty file as colliding'''
        self.assertRaises(RuntimeError, self.sub.check_for_collisions, self.empty_fn)

    def test_check_for_collision_empty(self):
        '''should identify full file as colliding'''
        self.assertRaises(RuntimeError, self.sub.check_for_collisions, self.full_fn)
        
    def test_check_for_collision_no(self):
        '''should identify destination as empty'''
        self.sub.check_for_collisions(self.no_fn)
        

def TestPipelineSteps(TestWithFiles):
    '''tests for file input/output for each pipeline step'''
    
    def test_make_otu_table(self):
        '''should read inputs and write otu table'''

        uc_lines = ["H  0  250  99.9  +  0  0  100M  seq0;counts=100  otu0", "H  0  250  99.9  +  0  0  100M  seq1;counts=10  otu0", "H  0  250  99.9  +  0  0  100M  seq5;counts=22  otu3"]
        index_lines = ['donor1  seq0  5', 'donor1 seq5 1', 'donor1 seq1 1', 'donor2  otu3  4']
        
        with open('tests/otus.uc', 'w') as f:
            for line in uc_lines:
                f.write(line + "\n")
                
        with open('tests/q.index', 'w') as f:
            for line in index_lines:
                f.write(line + "\n")
                
        subprocess.call(['python', 'uc2otus.py', 'tests/otus.uc', 'tests/q.index', '--output', 'tests/otu.count'])
        
        with open('tests/otus.count') as f:
            table_content = f.read()
            
        expected_table_data = [['OTU_ID', 'donor1', 'donor2'], ['otu1', '6', '0'], ['otu3', '1', '4']]
        expected_table = '\n'.join(['\t'.join(row) for row in table])
        
        self.assertEqual(table_content, expected_table)
        

if __name__ == '__main__':
    unittest.main(verbosity=2)
