### Optional preliminary steps
* Convert to new format
* Merge forward and reverse reads

### Middle steps
* Remove primers
* Quality filter (depends on removed primers fastq)
* Dereplicate (depends on quality filtered fasta)
* Make an index file (depends on dereplicated fasta and quality-filtered fasta)
* Make a sequence table (depends on quality filtered fasta)

### OTU calling
* Calling OTUs (depends on dereplicated fasta)
* Make an OTU table (depends on index file and uc files from OTU calling)