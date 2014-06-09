import unittest
from SmileTrain import uc2otus

class TestOTU(unittest.TestCase):
    '''tests for uc -> OTU'''
    
    def setUp(self):
        self.uc_lines = ["H       112678  253     99.6    +       0       0       520I53MD199M736I        seq0;counts=118114      2045", "H       108455  252     99.6    +       0       0       531I252M717I    seq1;counts=93322      38278"]
    
    def test_parse_uc_line(self):
        '''should parse one uc line into relevant fields'''
        line = self.uc_lines[0]
        hit, label, otu = uc2otus.parse_uc_line(line)
        self.assertEqual(hit, "H")
        self.assertEqual(label, "seq0;counts=118114")
        self.assertEqual(otu, "2045")
        
    def test_parse_uc_lines(self):
        '''should parse uc lines into dictionary'''
        d = uc2otus.parse_uc_lines(self.uc_lines)
        self.assertEqual(d, {'seq0': '2045', 'seq1': '38278'})
        
    def test_parse_bad_uc_lines(self):
        '''should throw exception when sequence ID repeated'''
        bad_uc_lines = ["H  0  250  99.9  +  0  0  100M  seq0;counts=100   otu1", "H  0  250  99.9  +  0  0  100M  seq0;counts=200  otu1"]
        self.assertRaises(RuntimeError, uc2otus.parse_uc_lines, bad_uc_lines)
        
    def test_parse_sample_lines(self):
        '''should read barcodes are first field'''
        bc_lines = ['animal4       TCCGTGCG', 'animal1       TCAAAGCT']
        self.assertEqual(uc2otus.parse_sample_lines(bc_lines), ['animal4', 'animal1'])
        
    def test_parse_index_lines(self):
        '''should split and recast index line'''
        line = 'donor1  seq5    100'
        self.assertEqual(uc2otus.parse_index_line(line), ['donor1', 'seq5', 100])
        
    def test_sparse_count_table(self):
        '''should make a sparse count table'''
        seq_otu = {'seqA': 'otuA', 'seqA2': 'otuA', 'seqC': 'otuC'}
        index_lines = ['donor1  seqA  5', 'donor1 seqA2  1', 'donor1 seqC 1', 'donor2  seqA  4']
        table = uc2otus.sparse_count_table(seq_otu, index_lines)
        self.assertEqual(table, {'donor1': {'otuA': 6, 'otuC': 1}, 'donor2': {'otuA': 4}})
        
    def test_otu_table(self):
        '''should make an OTU table'''
        table = {'donor1': {'otuA': 6, 'otuC': 1}, 'donor2': {'otuA': 4}}
        otus = ['otuA', 'otuC']
        samples = ['donor1', 'donor2']
        
        out_lines = [line for line in uc2otus.otu_table(table, otus=None, samples=samples)]
        
        self.assertEqual(out_lines[0], 'OTU_ID\tdonor1\tdonor2')
        self.assertIn('otuA\t6\t4', out_lines[1:])
        self.assertIn('otuC\t1\t0', out_lines[1:])
        
    def test_nomatch_first(self):
        '''should put no_match as first OTU ID if it's there'''
        table = {'donor1': {'otuA': 6, 'otuC': 1}, 'donor2': {'otuA': 4, 'no_match': 10}}
        out_lines = [line for line in uc2otus.otu_table(table)]
        first_id = out_lines[1].split()[0]
        
        self.assertEqual(first_id, 'no_match')