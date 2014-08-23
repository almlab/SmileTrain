## Was my sequencing run bad?

Some things to keep in mind here are:
* Were yours the only samples in the sequencing lane? If not, then you will expect that a smaller fraction of the total reads will map to your known barcodes.
* How many reads did you get? The forward and reverse ``.fastq`` files should be the same length. You can get its total length with ``wc -l`` and divide by 4 (but don't do this on the head node!); or, you can get the first 1000 entries (``head -4000 your.fastq``), check its size (``ls -lah``), and then extrapolate the number of reads in the big fastq.
  * New MISEQ runs should have about 25 million reads total. Does this include ones that go to phiX?
* What do your qualities look like? The BMC data includes a file ``...fastqc_report.html`` that shows how the quality of the read depends on nucleotide position.
* For MiSeq, you should be mostly in the green up to about 100 nucleotides.

## Santity checks

### After removing primers
* What fraction of your total reads were dropped?

### After barcoding
* How many index reads mapped to known barcodes?
* Were the most poorly-matched index reads the most lowly abundant?