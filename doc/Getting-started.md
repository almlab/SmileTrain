## Cloning
If you are going to be doing anything computationally intensive (which you probably are), you should put the SmileTrain code on your compute cluster (probably coyote).

Right now the Smile Train is not bundled as a package or anything, so you'll have to download the scripts and put them in the right place. I personally recommend `~/lib/SmileTrain`.

For example, if you are placing into `~/lib`, you should `cd ~/lib` and then `git clone https://github.com/almlab/SmileTrain.git`. (You can also get that https address from the main repository webpage.)

## Editing user.cfg
Before trying to run any of the scripts, you need to create a `user.cfg`. This file tells the SmileTrain scripts where to place temporary job submission scripts, where to look for the other scripts, etc.

A template is provided in the repository as `user.cfg.template`. You will definitely need to change the `username`, `tmp_directory`, `library`, and `bashrc` lines. (Make sure the `tmp_directory` folder exists!) The `queue` you pick will depend on your needs. (You can learn about the queues on your compute cluster with the obscure command `qmgr -c 'p s'` or the less informative `qstat -Q`.) You can point to my `usearch`, or you can download your own copy.

The `cluster` and `[Data]` lines are set up for use on coyote. If you are using a different cluster, you'll have to adjust those lines.

## Getting the right version of python
SmileTrain depends on some features of python that are specific to certain versions. You'll need python 2.7 (2.7.3 is the development version). You can see which version of python you are using by default by issuing `python --version`. If the version if not 2.7, you'll need to change it. On a cluster, this might mean manually calling `module load python/2.7.3` and/or adding that command to your `~/.bashrc`.

## Setting up your data
### From BMC-style raw Illumina data
You'll need forward and reverse reads (I'll call them `for.fastq` and `rev.fastq`) in [Illumina 1.3-1.7 format](http://drive5.com/usearch/manual/quality_score.html), a barcode file (I'll call it `barcode.txt`; it should have lines with sample name and barcode separated by a tab), and the forward and reverse primer (I'll call them AAA and TTT).

If you want to go all the way from your raw data to a reference-based OTU table using Greengenes, you'll just need to run `/path/to/SmileTrain/otu_caller.py -f for.fastq -r rev.fastq -p AAA -q TTT -b barcode.txt --all -n 10`.

The `--all` is a shortcut for `--check --split --convert --primers --merge --demultiplex --qfilter --dereplicate --index --ref_gg`. The `-n 10` means that the early steps (converting fastq format through quality filtering) will be performed in parallel on 10 nodes. The number of nodes you pick should be decided from a balance between job submission overhead and pure computational time needed.

### From Illumina 1.8 format
This is just the same as the above pipeline except that you don't need `--convert`, since you are already in the right file format. For example, if you want to map to Greengenes and make an OTU table, call `/path/to/SmileTrain/otu_caller.py -f for.fastq -r rev.fastq -p AAA -q TTT -b barcode.txt -n 10 --split --primers --merge --demultiplex --qfilter --dereplicate --index --ref_gg --otu_table`.

### From a fasta file
If you are starting with a QIIME fasta file, you should read [[How to process a filtered QIIME fasta]].

### Barcode mapping
The barcode mapping file is tab-delimited and has format `sample    barcode`, for example, `donor1_day1    ACGT`. Every sample-barcode combination goes on its own row. No headers por favor.

## It didn't work
Oh gosh, that's a different topic: [[Troubleshooting]]