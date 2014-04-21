intro to the pipeline
=======================================

Before doing anything, set up that config file.

then run the scripts in order

Input data
----------

Input data should be in FASTQ format. A FASTQ file has entries of four lines each:
	* the sequence ID (or "read ID") line, which is required to start with the symbol ``@``
	* the sequence line, which contains the sequence using the characters ACGT
	* a line that must begin with ``+`` and which can be empty of be a copy of the ID line
	* the quality line, whose alphabet can differ

New data from the BMC looks like this, so that's how I'm rolling right now.

::

	@MISEQ:1:1101:19196:1927#CTAGAATC/1
	TCGAAGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTATTGGGTGTAAAGGGTGCGTAGACGGCATGGCAAGTCTGATGTGAAAACCCGCGGCTCAACCCCGGCACTGCATTGCATCCTGCCAGCCTTGAGTGCCGGTGTGGCAAGTGGAATTCCTTGTGTACCGGTGAAATGCGTACATTTCCCGAGGAACTCCAGTTCCGAAGCCGGCTTCCTGCACGATCTCTGACGTTCT
	+MISEQ:1:1101:19196:1927#CTAGAATC/1
	P]]PP``cecP_`PPePPO`d``OecefePO``d^dN`baddd^efedbNdd_eQcfgddN]_QQQPNONNN^]NNP^dN][aOO[PNaaaPPaPaBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

This format includes the barcode read (or "index read" or "multiplex tag") on the ID line after a ``#``. This is (apparently) how the Illumina pipeline exports its data for Cassava versions 1.4 through 1.8, at which time the ID line will look like

::

	@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

where the barcode is now at the very end of the line.


Splitting FASTQs
----------------

to preserve memory


Removing primers
----------------

The input entry above was from a run with primer sequence

::

	ACGAAGTGCCAGCMGCCGCGGTAA

``remove_primers.py`` searches for this primer in the sequence line and then trims the sequence and quality lines to produce an output like

::

	@MISEQ:1:1101:19196:1927#CTAGAATC/1
	TACGTAGGGGGCAAGCGTTATCCGGATTTATTGGGTGTAAAGGGTGCGTAGACGGCATGGCAAGTCTGATGTGAAAACCCGCGGCTCAACCCCGGCACTGCATTGCATCCTGCCAGCCTTGAGTGCCGGTGTGGCAAGTGGAATTCCTTGTGTACCGGTGAAATGCGTACATTTCCCGAGGAACTCCAGTTCCGAAGCCGGCTTCCTGCACGATCTCTGACGTTCT
	+
	ecefePO``d^dN`baddd^efedbNdd_eQcfgddN]_QQQPNONNN^]NNP^dN][aOO[PNaaaPPaPaBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

The primer ended on the nucleotide just below the ``#`` on the previous line, so that's where the trimmed begins.


Demultiplexing
--------------