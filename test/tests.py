#!/usr/bin/env python

'''
Unit tests for the pipeline.

There is a child class TestWithFiles. This should be used as the parent for classes that
include tests that require reading and writing files. These files should be written into
the tests folder and then cleaned up.

Most of these tests should be refactored into separate scripts.
'''

from SmileTrain.test import fake_fh
import unittest, tempfile, subprocess, os, shutil, StringIO
from Bio import SeqIO, Seq

from SmileTrain import util, remove_primers, derep_fulllength, intersect, check_fastq_format, convert_fastq, map_barcodes, derep_fulllength, uc2otus, index

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
        
    def test_check_for_existence(self):
        '''should identify files as existing or not'''
        
        # the empty and full files should exist
        util.check_for_existence(self.empty_fn)
        util.check_for_existence(self.full_fn)
        
        # but the file we deleted shouldn't
        with self.assertRaises(RuntimeError):
            util.check_for_existence(self.no_fn)
            
    def test_check_for_collision_yes(self):
        '''should identify files as existing'''
        self.assertRaises(RuntimeError, util.check_for_collisions, self.empty_fn)
        self.assertRaises(RuntimeError, util.check_for_collisions, self.full_fn)
        
    def test_check_for_collision_no(self):
        '''should identify destination as empty'''
        util.check_for_collisions(self.no_fn)
        

class TestDereplicate(unittest.TestCase):
    '''tests for dereplication without samples'''
    
    def setUp(self):
        fasta = fake_fh(['>foo', 'AA', '>bar', 'CC', '>baz', 'AA', '>blag', 'AA', '>flog', 'TT', '>blob', 'TT'])
        minimum_counts = 2
        self.derep = derep_fulllength.Dereplicator(fasta, minimum_counts)
        
    def test_new_seq_id(self):
        '''should give increasing sequence ids'''
        it = self.derep.iter_seq_ids()
        self.assertEqual(it.next(), "seq0")
        self.assertEqual(it.next(), "seq1")
        self.assertEqual(it.next(), "seq2")
        
    def test_dereplicate_abundances(self):
        '''should properly count abundances'''
        self.assertEqual(self.derep.abundances, {'AA': 3, 'CC': 1, 'TT': 2})
        
    def test_dereplicate_seq_ids(self):
        '''should properly organize names'''
        self.assertEqual(self.derep.seq_ids, {'AA': 'seq0', 'CC': 'seq1', 'TT': 'seq2'})
        
    def test_sort_abundant_sequences(self):
        '''should properly short abundances'''
        self.assertEqual(self.derep.filtered_abundant_sequences, ['AA', 'TT'])
    
    def test_fasta_entries(self):
        '''should give abundance-sorted entries'''
        fe = self.derep.new_fasta_entries()
        record1 = fe.next()
        record2 = fe.next()
        self.assertEqual([record1.id, str(record1.seq), record2.id, str(record2.seq)], ['seq0;counts=3', 'AA', 'seq2;counts=2', 'TT'])
        

class TestIndex(unittest.TestCase):
    '''tests for index-writing script'''
        
    def test_parse_derep_fasta(self):
        '''should make a dictionary of fasta lines'''
        fasta = fake_fh(['>seq0;counts=10', 'AAA', '>seq4;counts=23', 'TTT'])
        self.assertEqual(index.parse_derep_fasta(fasta), {'AAA': 'seq0', 'TTT': 'seq4'})
        
    def test_sid_to_sample(self):
        '''should extract sample from fasta line'''
        self.assertEqual(index.sid_to_sample('sample=donor1;444'), 'donor1')
        
    def test_parse_full_fasta(self):
        seq_sid = {'AAA': 'seq0', 'TTT': 'seq4'}
        fasta = fake_fh(['>sample=donor1;1', 'AAA', '>sample=donor1;2', 'AAA', '>sample=donor1;3', 'TTT', '>sample=donorT;1', 'TTT'])
        abund = index.parse_full_fasta(fasta, seq_sid)
        self.assertEqual(abund, {('donor1', 'seq0'): 2, ('donor1', 'seq4'): 1, ('donorT', 'seq4'): 1})
        

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
