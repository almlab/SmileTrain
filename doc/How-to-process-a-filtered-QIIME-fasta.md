### The data
Imagine your collaborator has given you a file `lots_o_poo.fna`, which contains lines like

`>806rcbc0_0 M02171:6:000000000-A6FUV:1:1101:21715:4030 1:N:0:1 orig_bc=GGAGACAAGGGA new_bc=GGAGACAAGGGA bc_diffs=0
TACGGAGGGTCCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTACGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCAGGGCTCAACCTTGGAACTGCATT
`

This is a quality-filtered QIIME fasta file. The first line identifies the sample name (806rcbc0), the ID for this read in that sample (the 0 after the _), and a whole bunch of other information. The second line is the actual sequence.

### What the SmileTrain can do
Before boarding the Smile Train, you'll need to create a file `q.fst` that has usearch-style labels. For this example, it would look like `>sample=806rcbc0;0`.

There is a script `SmileTrain/util/qiime_to_st_labels.py` that can do this for you! Make sure you are not on the head node (maybe `qsub -I`) and then run

`/path/to/SmileTrain/util/qiime_to_st_labels.py /path/to/lots_o_poo.fna > /path/to/q.fst`

Take a look to see that `q.fst` looks the way you expect. Get back on the head node, and then you can get greengenes OTUs and some OTU tables by running

`/path/to/SmileTrain/otu_caller.py --dereplicate --index --ref_gg --otu_table`

Et voilà.