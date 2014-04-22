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

Primer removal and quality filtering are both embarrassingly parallel computing tasks: every FASTQ entry can be processed independent of all the others. The pipeline might split up the FASTQ files to utilize all available CPUs for this part of the job.

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

See also: :doc:`primers`

Merging reads
-------------

Forward and reverse reads that came from the same molecule should have the same sequence ID. If they are long enough to overlap, they can be *intersected* or *merged*. In some cases, the merged read allows improved confidence of nucleotide values in the middle of the read where qualities are typically most poor.

The pipeline currently supports merging with usearch. It should aim to also use Sonia's SHE-RA.

Demultiplexing
--------------

Every forward and reverse read has a barcode associated with it. If a read doesn't match a known barcode, that read should be thrown out. If a barcode matches a known barcode, the read is assigned to that sample.

For example, if there is a line in the barcode mapping file

::

	HSD_group2_day15_animal7	CTAGAATC

then we can tell that the above read came from that sample, so after demultiplexing that FASTQ entry will read

::

	@sample=HSD_group2_day15_animal7;234
	TACGTAGGGGGCAAGCGTTATCCGGATTTATTGGGTGTAAAGGGTGCGTAGACGGCATGGCAAGTCTGATGTGAAAACCCGCGGCTCAACCCCGGCACTGCATTGCATCCTGCCAGCCTTGAGTGCCGGTGTGGCAAGTGGAATTCCTTGTGTACCGGTGAAATGCGTACATTTCCCGAGGAACTCCAGTTCCGAAGCCGGCTTCCTGCACGATCTCTGACGTTCT
	+
	ecefePO``d^dN`baddd^efedbNdd_eQcfgddN]_QQQPNONNN^]NNP^dN][aOO[PNaaaPPaPaBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

where ``234`` indicates that this was the 234th read associated with that sample. Only the first line has changed.

Quality filtering
-----------------

Every nucleotide in a read has an associated quality score. If we have low confidence in a nucleotide, we should throw it out. Reads that have poor quality overall should also be thrown out.

The pipeline currently supports quality filtering using usearch.

Dereplication
-------------

OTU calling
-----------

Reference-based
~~~~~~~~~~~~~~~

De novo
~~~~~~~

Making OTU tables
-----------------

