import unittest, tempfile, subprocess, os

import util, derep_fulllength

class TestFastaUtilities(unittest.TestCase):
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

		os.rmdir('tests')


if __name__ == '__main__':
    unittest.main(verbosity=2)