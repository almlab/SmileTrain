import unittest
from SmileTrain import seq_table
from SmileTrain.test import FakeFile

class TestFastaToTableAndAbundance(unittest.TestCase):
    def test_correct(self):
        fasta = FakeFile(['>sample=donor1;1', 'AAA', '>sample=donor1;2', 'AAA', '>sample=donor2;1', 'AAA', '>sample=donor3;1', 'TTT'])
        table = {'AAA': {'donor1': 2, 'donor2': 1}, 'TTT': {'donor3': 1}}
        abund = {'AAA': 3, 'TTT': 1}
        self.assertEqual(seq_table.fasta_to_table_and_abund(fasta), (table, abund))
    
class TestTableToSamples(unittest.TestCase):
    def test_correct(self):
        table = {"AAA": {"donor1": 5, "donor2": 10}, "CCC": {"donor1": 1, "donor5": 4}, "TTT": {"donor4": 1}}
        samples = ['donor1', 'donor2', 'donor4', 'donor5']
        self.assertEqual(seq_table.table_to_samples(table), samples)
    
class TestTableLines(unittest.TestCase):
    def test_correct(self):
        table = {"AAA": {"donor1": 5, "donor2": 10}, "CCC": {"donor1": 1, "donor5": 4}, "TTT": {"donor4": 1}}
        lines = ['sequence\tdonor1\tdonor2\tdonor4\donor5', 'AAA\t5\t10\t0\t0', 'CCC\t1\t0\t0\t4', 'TTT\t0\t0\t1\t0']