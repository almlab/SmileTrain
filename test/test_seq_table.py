import unittest
from SmileTrain import index_to_seq_table

class TestTrimmedFastaEntries(unittest.TestCase):
    def test_correct(self):
        lines = ['>seq1;counts=500', 'AAA', '>seq2;counts=123', 'GGG']
        self.assertEqual(index_to_seq_table.trimmed_fasta_entries(lines), [['seq1', 'AAA'], ['seq2', 'GGG']])
     
   
class TestSparseSeqCountTable(unittest.TestCase):    
    def test_correct(self):
        lines = ['COSM2T1C	seq861547	1', 'COSM3T2oilsCsP	seq392727	1', 'GAB14NB	seq1300749	1', 'COSM2T1C	seq795681	3']
        table = {'COSM2T1C': {'seq861547': 1, 'seq795681': 3}, 'COSM3T2oilsCsP': {'seq392727': 1}, 'GAB14NB': {'seq1300749': 1}}
        self.assertEqual(index_to_seq_table.sparse_seq_count_table(lines), table)


class TestSeqTable(unittest.TestCase):    
    def test_correct(self):
        table = {'COSM2T1C': {'seq861547': 1, 'seq795681': 3}, 'COSM3T2oilsCsP': {'seq392727': 1}, 'GAB14NB': {'seq1300749': 1}}
        lab_seq = [['seq1300749', 'TTT'], ['seq861547', 'AAA'], ['seq795681', 'CCC'], ['seq392727', 'GGG']]
        lines = ['SEQUENCE\tCOSM2T1C\tCOSM3T2oilsCsP\tGAB14NB', 'TTT\t0\t0\t1', 'AAA\t1\t0\t0', 'CCC\t3\t0\t0', 'GGG\t0\t1\t0']
        out = list(index_to_seq_table.seq_table(table, lab_seq))
        self.assertEqual(lines, out)