#!/usr/bin/env python

import unittest, tempfile, subprocess, os, shutil

import util, remove_primers, derep_fulllength, intersect, check_fastq_format

class TestWithFiles(unittest.TestCase):
    
    def tearDown(self):
        if os.path.isdir('tests'):
            shutil.rmtree('tests')


class TestFastaUtilities(TestWithFiles):
    '''tests for scripts and functions that do simple fasta manipulations'''

    def setUp(self):
        self.good_fasta_content = ">foo\nAAA\nAAA\n>bar\nCCC\n>baz\nGGG"
        self.bad_fasta_content = ">foo\nAAA\nAAA\n>bar\n>baz\nGGG"

    def test_fasta_entries_with_good_content(self):
        '''fasta_entries should split the content as expected'''

        fasta_lines = self.good_fasta_content.split('\n')
        i = util.fasta_entries(fasta_lines)

        self.assertEqual(i.next(), ["foo", "AAAAAA"])
        self.assertEqual(i.next(), ["bar", "CCC"])
        self.assertEqual(i.next(), ["baz", "GGG"])

    def test_fasta_entries_with_bad_content(self):
        '''fasta_entries should raise AssertionError when it gets two > lines back to back'''

        fasta_lines = self.bad_fasta_content.split('\n')
        i = util.fasta_entries(fasta_lines)

        i.next()
        with self.assertRaises(AssertionError):
            i.next()

    def test_split_fasta(self):
        '''split_fasta.py should split content as expected'''

        os.mkdir('tests')
        fasta_fh, fasta_fn = tempfile.mkstemp(suffix='.fasta', dir='tests')
        os.write(fasta_fh, self.good_fasta_content)
        os.close(fasta_fh)

        subprocess.call(['python', 'split_fasta.py', fasta_fn, '3'])
        with open(fasta_fn + '.0', 'r') as f:
            fasta_out0 = f.read()

        self.assertEqual(fasta_out0, ">foo\nAAAAAA\n")

        os.remove(fasta_fn)
        for i in range(3):
            os.remove("%s.%d" %(fasta_fn, i))


class TestFastqUtilities(TestWithFiles):
    '''tests for functions and scripts that do fastq manipulations'''
    def setUp(self):
        self.good_fastq_content = "@foo\nAAA\n+foo\n!!!\n@bar\nCCC\n+bar\n###"
        
    def test_split_fastq(self):
        '''split_fastq.py should split and trim content as expected'''
        
        os.mkdir('tests')
        fastq_fh, fastq_fn = tempfile.mkstemp(suffix='.fastq', dir='tests')
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
    
            
class TestFileChecks(TestWithFiles):
    '''tests for utilities that make queries about files'''
    
    def setUp(self):
        os.mkdir('tests')
        
        empty_fh, self.empty_fn = tempfile.mkstemp(dir='tests')
        
        full_fh, self.full_fn = tempfile.mkstemp(dir='tests')
        os.write(full_fh, "hello world")
        os.close(full_fh)
        
        no_fh, self.no_fn = tempfile.mkstemp(dir='tests')
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
            
    def test_check_for_collision(self):
        '''should identify files and not existing or so'''


class TestRemovePrimers(unittest.TestCase):
    '''tests for the remove primers step'''

    def setUp(self):
        self.fastq = "@lolapolooza\nTAAAACATCATCATCAT\n+whatever\nabcdefghijklmnopq\n"
        self.fastq_lines = self.fastq.split()
        self.primer = "AAAA"
        self.max_primer_diffs = 1

        self.primer_remover = remove_primers.PrimerRemover(self.fastq_lines, self.primer, self.max_primer_diffs)

    def test_correct_output(self):
        '''the primer remover should trim the match as expected'''
        self.assertEqual(self.primer_remover.next(), "@lolapolooza\nCATCATCATCAT\n+\nfghijklmnopq")
        self.assertEqual(self.primer_remover.n_successes, 1)
        
        
class TestIntersect(TestWithFiles):
    '''tests for merging reads'''
    
    def test_usearch_merge(self):
        '''usearch should merge a perfectly matched pair of sequences'''
        
        # swo> I choose 300 random ACGTs, then took the first ~200 and last ~200, reverse
        # complemented the last, gave them all the same quality code, and made up a name.
        self.forward_fastq = "@foo/1\nGTTTTCTTCGCTTTATGGTGGTGGTAAAAGTGCTTCGATCTGCTAGATATCCCTCAGGAAAGTTTATGCCCGTGTCCGTTTGTTTGGGTAGATCTCTCACCCTTGGAATTCCAAGCGTTCAGGTATCCCACAATCGCTTCGATGACTCCGCCTCCTTATTATATACTTCGCCGATACGCAGCGCATGAAGAGTCATCGGGA\n+\n#########################################################################################################################################################################################################\n"
        self.reverse_fastq = "@foo/2\nCGATATCCGTGGCTTAAGCTATATGCGATTTTGCAGAGCAGTCAAGGTCTCCCTGGGTAGATTAAAGGGCGAGCTCACGAAGAGATTACTACTCAACCCTCCCGATGACTCTTCATGCGCTGCGTATCGGCGAAGTATATAATAAGGAGGCGGAGTCATCGAAGCGATTGTGGGATACCTGAACGCTTGGAATTCCAAGG\n+\n########################################################################################################################################################################################################\n"
        
        # set up the tests directory
        os.mkdir('tests')
        
        fasta_fh, fasta_fn = tempfile.mkstemp(suffix='.fasta', dir='tests')
        
        with open('tests/f.fastq', 'w') as f:
            f.write(self.forward_fastq)
            
        with open('tests/r.fastq', 'w') as f:
            f.write(self.reverse_fastq)

        # run usearch. usearch should output something that has "1  Exact overlaps" in it.
        usearch_output = subprocess.check_output(['usearch', '-fastq_mergepairs', 'tests/f.fastq', '-reverse', 'tests/r.fastq', '-fastqout', 'tests/out.fastq'], stderr=subprocess.STDOUT)
        self.assertRegexpMatches(usearch_output, '1\s+Exact overlaps')
        
        # check the output too
        with open('tests/out.fastq', 'r') as f:
            out_fastq = f.read()
            
        self.assertEqual(out_fastq, "@foo/1\nGTTTTCTTCGCTTTATGGTGGTGGTAAAAGTGCTTCGATCTGCTAGATATCCCTCAGGAAAGTTTATGCCCGTGTCCGTTTGTTTGGGTAGATCTCTCACCCTTGGAATTCCAAGCGTTCAGGTATCCCACAATCGCTTCGATGACTCCGCCTCCTTATTATATACTTCGCCGATACGCAGCGCATGAAGAGTCATCGGGAGGGTTGAGTAGTAATCTCTTCGTGAGCTCGCCCTTTAATCTACCCAGGGAGACCTTGACTGCTCTGCAAAATCGCATATAGCTTAAGCCACGGATATCG\n+\n####################################################################################################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$###################################################################################################\n")


if __name__ == '__main__':
    unittest.main(verbosity=2)