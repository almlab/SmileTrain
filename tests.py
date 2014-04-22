#!/usr/bin/env python

import unittest, tempfile, subprocess, os, shutil

import util, remove_primers, derep_fulllength, intersect

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
    
    def test_fastq_id_parsing(self):
        '''fastq_at_line_to_id should trim the starting @ and trailing /1 or /2'''
        
        self.assertEqual(intersect.fastq_at_line_to_id("@lolapolooza:1234#ACGT/1"), "lolapolooza:1234#ACGT")
        self.assertEqual(intersect.fastq_at_line_to_id("@lolapolooza:1234#ACGT/2"), "lolapolooza:1234#ACGT")
        
    def test_fastq_ids(self):
        '''fastq_ids should get the right ids from some fastq entries'''
        
        fastq1 = "@lolapolooza:1234#ACGT/1\nTAAAACATCATCATCAT\n+whatever\nabcdefghijklmnopq\n"
        fastq2 = "@lolapolooza:7890#TGCA/1\nGAATACTACGGGAGAGAAA\n+whatever\nabcdefghijklmnopqrs\n"
        fastq_lines = (fastq1 + fastq2).split('\n')
        
        ids = intersect.fastq_ids(fastq_lines)
        self.assertEqual(ids, ["lolapolooza:1234#ACGT", "lolapolooza:7890#TGCA"])
    
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