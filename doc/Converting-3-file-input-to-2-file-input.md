## What are 2- and 3-file inputs?
Illumina sequencing centers other than the BMC give you three sequence files: a forward fastq, a reverse fastq, and an index fastq. The Smile Train is currently set up to leave Two-File Station, which means that you'll need to put the index reads into the headers of the forward and reverse reads.

A three-file format has a forward fastq with entries like
```
@M02171:14:000000000-A6CAE:1:1101:15953:1560 1:N:0:101
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTTGTT
+
BBBBBBBBBBBBGGGGGGGGGGHGGGGGHHHHFFHHGGGGCDGG0EEGGGGGHHHGH
```
and reverse fastq like
```
@M02171:14:000000000-A6CAE:1:1101:15953:1560 2:N:0:101
CCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGCA
+
ABCCCFFFFFFFGGGGGG2FFGGGGGGGHHH5FH5DFHHHHHHHFHFHHHG
```
and an index fastq like
```
@M02171:14:000000000-A6CAE:1:1101:15953:1560 1:N:0:101
ACCACATACATC
+
FFFFFFFFFFFF
```
It is the number (`1` or `2`) in the front of the last few fields in the `@` line that specifies the read direction (i.e., `1:N:0:101` for forward or index, `2:N:0:101` for reverse).

In a two-file format, the index read is put right into the `@` line. The corresponding two-file entries for the three-file entries above would be
```
@M02171:14:000000000-A6CAE:1:1101:15953:1560#ACCACATACATC/1
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTTGTT
+
BBBBBBBBBBBBGGGGGGGGGGHGGGGGHHHHFFHHGGGGCDGG0EEGGGGGHHHGH
```
and
```
@M02171:14:000000000-A6CAE:1:1101:15953:1560#ACCACATACATC/2
CCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGCA
+
ABCCCFFFFFFFGGGGGG2FFGGGGGGGHHH5FH5DFHHHHHHHFHFHHHG
```
Note that the `@` line has no whitespace, the barcode read comes after the `#` and before the `/`, and either `1` or `2` follows the `/` to indicate read direction.

## What scripts to run
There is a utility script that will do this. To convert `for.fastq`, `rev.fastq`, and `ind.fastq` into new fastqs `new_for.fastq` and `new_rev.fastq`, issue

`/path/to/SmileTrain/tools/convert_3file_to_2file.py for.fastq rev.fastq ind.fastq new_for.fastq new_rev.fastq`

This can be a little slow, so be sure to _not_ run it on the head node.

## What if there is some complaint about unequal quality and sequence lengths?
I've seen some cases where the all the index read entries have a quality score line that has a different length from the sequence line. If you think this is just an artifact of the Illumina pipeline and you actually trust those index reads, you can run `fix_index_fastq.py` to either pad or trim the quality score line as necessary.

`/path/to/SmileTrain/tools/fix_index_fastq.py old_index.fastq > new_index.fastq`