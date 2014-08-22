import unittest
from SmileTrain import util_primer
from SmileTrain.test import fake_fh

class TestRemovePrimers(unittest.TestCase):
    '''tests for the remove primers utility'''

    def setUp(self):
        self.fastq = fake_fh('''@lolapolooza\nTAAAACATCATCATCAT\n+lolapolooza\n"#$%&'()*+,-./012\n''')
        self.primer = "AAAA"
        self.max_primer_diffs = 1

        self.primer_remover = util_primer.PrimerRemover(self.fastq, self.primer, self.max_primer_diffs)

    def test_correct_output(self):
        '''the primer remover should trim the match as expected'''
        record = self.primer_remover.next()
        self.assertEqual(record.id, 'lolapolooza')
        self.assertEqual(str(record.seq), 'CATCATCATCAT')
        self.assertEqual(record.letter_annotations['phred_quality'], [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])
        self.assertEqual(self.primer_remover.n_successes, 1)
     
        
class TestMistmatches(unittest.TestCase):
    def test_correct(self):
        self.assertEqual(util_primer.mismatches('TCAAAAGATGATGATGAT', 'AAAA', 15), (2, 0))