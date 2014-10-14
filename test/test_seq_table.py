import unittest
from SmileTrain import seq_table
from SmileTrain.test import fake_fh

class TestFastaToTableAndAbundance(unittest.TestCase):
    def test_correct(self):
        fasta = fake_fh(['>sample=donor1;1', 'AAA', '>sample=donor1;2', 'AAA', '>sample=donor2;1', 'AAA', '>sample=donor3;1', 'TTT'])
        names = {'AAA': 'seqA', 'TTT': 'seqT'}
        table = {'seqA': {'donor1': 2, 'donor2': 1}, 'seqT': {'donor3': 1}}
        abund = {'seqA': 3, 'seqT': 1}
        table2, abund2 = seq_table.SeqTableWriter.fasta_to_abund(fasta, names)
        self.assertEqual(table, table2)
        self.assertEqual(abund, abund2)
    
class TestTableToSamples(unittest.TestCase):
    def test_correct(self):
        table = {"AAA": {"donor1": 5, "donor2": 10}, "CCC": {"donor1": 1, "donor5": 4}, "TTT": {"donor4": 1}}
        samples = ['donor1', 'donor2', 'donor4', 'donor5']
        self.assertEqual(seq_table.SeqTableWriter.table_to_samples(table), samples)
    
class TestTableLines(unittest.TestCase):
    def test_correct(self):
        table = {"seq0": {"donor1": 5, "donor2": 10}, "seq1": {"donor1": 1, "donor5": 4}, "seq2": {"donor4": 1}}
        abund = {'seq0': 15, 'seq1': 5, 'seq2': 1}
        samples = ['donor1', 'donor2', 'donor4', 'donor5']
        cont = "\n".join(['sequence_id\tdonor1\tdonor2\tdonor4\tdonor5', 'seq0\t5\t10\t0\t0', 'seq1\t1\t0\t0\t4', 'seq2\t0\t0\t1\t0']) + "\n"
        fh = fake_fh()
        seq_table.SeqTableWriter.write_table(table, abund, samples, 0, fh)
        fh.seek(0)
        cont2 = fh.read()
        self.assertEqual(cont, cont2)