import unittest, tempfile, subprocess, os, shutil
from SmileTrain.tools import split_sequence_barcodes
from SmileTrain.test import fake_fh


class TestBestPlateBarcodeMatch(unittest.TestCase):
    def test_correct(self):
        '''should pick out a plate barcode from a sequence'''
        barcode = "YRYR"
        seq = "AACACATTTTTTTTTTTTTTTTTTT"
        bests = split_sequence_barcodes.best_plate_barcode_match([barcode], seq)
        self.assertEqual(bests, ('YRYR', 3, 0))


if __name__ == '__main__':
    unittest.main(verbosity=2)
