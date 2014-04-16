import unittest

class TestPipeline(unittest.TestCase):
	def setUp(self):
		self.fasta = "asdf"

	# make sure that the primer remover is getting rid of the right thing
	def test_remove_primer(self):