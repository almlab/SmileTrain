#!/usr/bin/env python

'''
Unit tests for the pipeline.

There is a child class TestWithFiles. This should be used as the parent for classes that
include tests that require reading and writing files. These files should be written into
the tests folder and then cleaned up.

Most of these tests should be refactored into separate scripts.
'''

import unittest, tempfile, subprocess, os, shutil

from SmileTrain import util, remove_primers, derep_fulllength, intersect, check_fastq_format, convert_fastq, map_barcodes, derep_fulllength, uc2otus, index

tmp_dir = 'tmp'

class TestWithFiles(unittest.TestCase):
    '''tests that need to read and write files'''
    
    def setUp(self):
        os.mkdir(tmp_dir)
    
    def tearDown(self):
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)


class TestFastaUtilities(TestWithFiles):
    '''tests for scripts and functions that do simple fasta manipulations'''

    def setUp(self):
        self.good_fasta_content = ">foo\nAAA\nAAA\n>bar\nCCC\n>baz\nTTT"
        self.bad_fasta_content = ">foo\nAAA\nAAA\n>bar\n>baz\nTTT"

    def test_fasta_entries_with_good_content(self):
        '''fasta_entries should split the content as expected'''

        fasta_lines = self.good_fasta_content.split('\n')
        i = util.fasta_entries(fasta_lines)

        self.assertEqual(i.next(), ["foo", "AAAAAA"])
        self.assertEqual(i.next(), ["bar", "CCC"])
        self.assertEqual(i.next(), ["baz", "TTT"])

    def test_fasta_entries_with_bad_content(self):
        '''fasta_entries should raise AssertionError when it gets two > lines back to back'''

        fasta_lines = self.bad_fasta_content.split('\n')
        i = util.fasta_entries(fasta_lines)

        i.next()
        with self.assertRaises(AssertionError):
            i.next()

def TestSplitFasta(TestFastaUtilities):
    def setUp(self):
        fasta_content = ">foo\nAAA\nAAA\n>bar\nCCC\n>baz\nTTT"
        fasta_fh, self.fasta_fn = tempfile.mkstemp(suffix='.fasta', dir=tmp_dir)
        os.write(fasta_fh, self.good_fasta_content)
        os.close(fasta_fh)
        
    def test_split_fasta(self):
        '''split_fasta.py should split content as expected'''

        subprocess.call(['python', 'split_fasta.py', self.fasta_fn, '3'])
        with open(self.fasta_fn + '.0', 'r') as f:
            self.assertEqual(f.read(), ">foo\nAAAAAA\n")
            
        with open(self.fasta_fn + '.1') as f:
            self.assertEqual(f.read(), ">bar\nCCC\n")
            
        with open(self.fasta_fn + '.2') as f:
            self.assertEqual(f.read(), ">baz\nTTT\n")

        os.remove(self.fasta_fn)
        for i in range(3):
            os.remove("%s.%d" %(self.fasta_fn, i))
            
    def test_split_fasta_by_hash(self):
        '''split_fasta should split sequences as expected'''
        
        subprocess.call(['python', 'split_fasta.py', '--hash', self.fasta_fn, '2'])
        
        with open(self.fasta_fn + '.0') as f:
            self.assertEqual(f.read(), ">foo\nAAAAAA\n>bar\nCCC\n")
            
        with open(self.fasta_fn + '.1') as f:
            self.assertEqual(f.read(), ">baz\nTTT\n")


class TestFastqUtilities(TestWithFiles):
    '''tests for functions and scripts that do fastq manipulations'''
    def setUp(self):
        os.mkdir(tmp_dir)
        self.good_fastq_content = "@foo\nAAA\n+foo\n!!!\n@bar\nCCC\n+bar\n###"
        self.fastq13 = """@lolapolooza:1234#ACGT/1\nAATTAAGTCAAATTTGGCCTGGCCCAGTGTCCAATGTTGT\n+\nABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"""
        self.fastq18 = """@lolapolooza:1234#ACGT/1\nAATTAAGTCAAATTTGGCCTGGCCCAGTGTCCAATGTTGT\n+\n"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHJ"""
        
    def test_fastq_iterator(self):
        it = util.fastq_iterator(self.good_fastq_content.split())
        self.assertEqual(it.next(), ['@foo', 'AAA', '!!!'])
        self.assertEqual(it.next(), ['@bar', 'CCC', '###'])
        
    def test_split_fastq(self):
        '''split_fastq.py should split and trim content as expected'''
        
        fastq_fh, fastq_fn = tempfile.mkstemp(suffix='.fastq', dir=tmp_dir)
        os.write(fastq_fh, self.good_fastq_content)
        os.close(fastq_fh)
        
        subprocess.call(['python', 'split_fastq.py', fastq_fn, '2'])
        
        with open("%s.0" % fastq_fn) as f:
            fastq_out0 = f.read()
        # the plus line has been expunged
        self.assertEqual(fastq_out0, "@foo\nAAA\n+\n!!!\n")
        
        with open("%s.1" % fastq_fn) as f:
            fastq_out1 = f.read()    
        self.assertEqual(fastq_out1, "@bar\nCCC\n+\n###\n")
        
        os.remove(fastq_fn)
        for i in range(2):
            os.remove("%s.%d" % (fastq_fn, i))
            
    def test_fastq_id_parsing(self):
        '''fastq_at_line_to_id should trim the starting @ and trailing /1 or /2'''
        
        self.assertEqual(util.fastq_at_line_to_id("@lolapolooza:1234#ACGT/1"), "lolapolooza:1234#ACGT")
        self.assertEqual(util.fastq_at_line_to_id("@lolapolooza:1234#ACGT/2"), "lolapolooza:1234#ACGT")
        
    def test_fastq_ids(self):
        '''fastq_ids should get the right ids from some fastq entries'''
        
        fastq1 = "@lolapolooza:1234#ACGT/1\nTAAAACATCATCATCAT\n+whatever\nabcdefghijklmnopq\n"
        fastq2 = "@lolapolooza:7890#TGCA/1\nGAATACTACGGGAGAGAAA\n+whatever\nabcdefghijklmnopqrs\n"
        fastq_lines = (fastq1 + fastq2).split('\n')
        
        ids = intersect.fastq_ids(fastq_lines)
        self.assertEqual(ids, ["lolapolooza:1234#ACGT", "lolapolooza:7890#TGCA"])
        
    def test_fastq_id_match(self):
        fastq1 = "@lolapolooza:1234#ACGT/1\nTAAAACATCATCATCAT\n+whatever\nabcdefghijklmnopq\n"
        fastq2 = "@lolapolooza:7890#TGCA/1\nGAATACTACGGGAGAGAAA\n+whatever\nabcdefghijklmnopqrs\n"
        fastq_lines = (fastq1 + fastq2).split('\n')
        
        out = intersect.fastq_entries_with_matching_ids(fastq_lines, ['lolapolooza:7890#TGCA']).next()
        self.assertEqual(out, '@lolapolooza:7890#TGCA/1\nGAATACTACGGGAGAGAAA\n+\nabcdefghijklmnopqrs')
        
    def test_format_identification(self):
        '''should discriminate Illumina 1.3-1.7 against 1.8 against trash'''
        fastq13_quality = "AJ[h"
        fastq18_quality = "AJ';"
        ambiguous_fastq = "ABCJ"
        mixed_fastq = "AJ[h';"
        bad_fastq = "AJ{|"
        
        self.assertEqual(check_fastq_format.quality_line_format(fastq13_quality), 'illumina13')
        self.assertEqual(check_fastq_format.quality_line_format(fastq18_quality), 'illumina18')
        self.assertEqual(check_fastq_format.quality_line_format(ambiguous_fastq), 'ambiguous')
        
        with self.assertRaises(RuntimeError):
            check_fastq_format.quality_line_format(mixed_fastq)
            
        with self.assertRaises(RuntimeError):
            check_fastq_format.quality_line_format(bad_fastq)
            
    def test_format_conversion(self):
        '''should convert Illumina 1.4-1.7 to our mixed format'''
        
        entry13 = util.fastq_iterator(self.fastq13.split("\n")).next()
        entry18 = util.fastq_iterator(self.fastq18.split("\n")).next()
        
        self.assertEqual(convert_fastq.convert_entry(entry13), entry18)
        
    def test_format_check_right(self):
        '''should identify written file as correct format'''
        with open('%s/good.fastq' % tmp_dir, 'w') as f:
            f.write(self.fastq13)
        
        check_fastq_format.check_illumina_format(['%s/good.fastq' % tmp_dir], 'illumina13')
        
    def test_format_check_wrong(self):
        '''should identify written file with incorrect format'''
        with open('%s/bad.fastq' % tmp_dir, 'w') as f:
            f.write(self.fastq18)
        
        self.assertRaises(RuntimeError, check_fastq_format.check_illumina_format, ['%s/bad.fastq' % tmp_dir], 'illumina13')
    
            
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
        
        
class TestIntersect(TestWithFiles):
    '''tests for merging reads'''
    
    def test_usearch_merge(self):
        '''usearch should merge a perfectly matched pair of sequences'''
        
        # swo> I choose 300 random ACGTs, then took the first ~200 and last ~200, reverse
        # complemented the last, gave them all the same quality code, and made up a name.
        self.forward_fastq = "@foo/1\nGTTTTCTTCGCTTTATGGTGGTGGTAAAAGTGCTTCGATCTGCTAGATATCCCTCAGGAAAGTTTATGCCCGTGTCCGTTTGTTTGGGTAGATCTCTCACCCTTGGAATTCCAAGCGTTCAGGTATCCCACAATCGCTTCGATGACTCCGCCTCCTTATTATATACTTCGCCGATACGCAGCGCATGAAGAGTCATCGGGA\n+\n#########################################################################################################################################################################################################\n"
        self.reverse_fastq = "@foo/2\nCGATATCCGTGGCTTAAGCTATATGCGATTTTGCAGAGCAGTCAAGGTCTCCCTGGGTAGATTAAAGGGCGAGCTCACGAAGAGATTACTACTCAACCCTCCCGATGACTCTTCATGCGCTGCGTATCGGCGAAGTATATAATAAGGAGGCGGAGTCATCGAAGCGATTGTGGGATACCTGAACGCTTGGAATTCCAAGG\n+\n########################################################################################################################################################################################################\n"
        
        fasta_fh, fasta_fn = tempfile.mkstemp(suffix='.fasta', dir=tmp_dir)
        
        with open('%s/f.fastq' % tmp_dir, 'w') as f:
            f.write(self.forward_fastq)
            
        with open('%s/r.fastq' % tmp_dir, 'w') as f:
            f.write(self.reverse_fastq)

        # run usearch. usearch should output something that has "1  Exact overlaps" in it.
        usearch_output = subprocess.check_output(['usearch', '-fastq_mergepairs', '%s/f.fastq' % tmp_dir, '-reverse', '%s/r.fastq' % tmp_dir, '-fastqout', '%s/out.fastq' % tmp_dir], stderr=subprocess.STDOUT)
        self.assertRegexpMatches(usearch_output, '1\s+Exact overlaps')
        
        # check the output too
        with open('%s/out.fastq' % tmp_dir, 'r') as f:
            out_fastq = f.read()
            
        self.assertEqual(out_fastq, "@foo/1\nGTTTTCTTCGCTTTATGGTGGTGGTAAAAGTGCTTCGATCTGCTAGATATCCCTCAGGAAAGTTTATGCCCGTGTCCGTTTGTTTGGGTAGATCTCTCACCCTTGGAATTCCAAGCGTTCAGGTATCCCACAATCGCTTCGATGACTCCGCCTCCTTATTATATACTTCGCCGATACGCAGCGCATGAAGAGTCATCGGGAGGGTTGAGTAGTAATCTCTTCGTGAGCTCGCCCTTTAATCTACCCAGGGAGACCTTGACTGCTCTGCAAAATCGCATATAGCTTAAGCCACGGATATCG\n+\n####################################################################################################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$###################################################################################################\n")
        

class TestDereplicate(unittest.TestCase):
    '''tests for dereplication without samples'''
    
    def setUp(self):
        entries = [['foo', 'AA'], ['bar', 'CC'], ['baz', 'AA'], ['blag', 'AA'], ['flog', 'TT'], ['blob', 'TT']]
        minimum_counts = 2
        self.derep = derep_fulllength.Dereplicator(entries, minimum_counts)
        
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
        self.assertEqual(fe.next(), '>seq0;counts=3\nAA')
        self.assertEqual(fe.next(), '>seq2;counts=2\nTT')
        

class TestIndex(unittest.TestCase):
    '''tests for index-writing script'''
    
    def test_parse_seq_id(self):
        '''should get sequence ID from first fasta line'''
        self.assertEqual(index.parse_seq_sid('seq44;counts=12'), 'seq44')
        
    def test_parse_derep_fasta(self):
        '''should make a dictionary of fasta lines'''
        fasta_lines = ['>seq0;counts=10', 'AAA', '>seq4;counts=23', 'TTT']
        self.assertEqual(index.parse_derep_fasta(fasta_lines), {'AAA': 'seq0', 'TTT': 'seq4'})
        
    def test_sid_to_sample(self):
        '''should extract sample from fasta line'''
        self.assertEqual(index.sid_to_sample('sample=donor1;444'), 'donor1')
        
    def test_parse_full_fasta(self):
        seq_sid = {'AAA': 'seq0', 'TTT': 'seq4'}
        fasta_lines = ['>sample=donor1;1', 'AAA', '>sample=donor1;2', 'AAA', '>sample=donor1;3', 'TTT', '>sample=donorT;1', 'TTT']
        abund = index.parse_full_fasta(fasta_lines, seq_sid)
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
