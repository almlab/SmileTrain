'''
Demultiplex reads by mapping the barcode read to the sample names from a barcode mapping
file.

Given a tab-separated barcode mapping file like
    ACGT    sample1 

input like
    @OURSEQ:lolapalooza1234#ACGT/1
    AACCGGTT
    +
    abcdefgh

becomes output like
    @sample=sample1;1
    AACCGGTT
    +whatever
    abcdefgh
'''

import usearch_python.primer
import sys, argparse, string, itertools, re

# swo> this code should probably be moved into a separate module
from remove_primers import mismatches

rctab = string.maketrans('ACGTacgt','TGCAtgca')
def reverse_complement(x):
    return x[::-1].translate(rctab)

def barcode_file_to_dictionary(barcode_lines):
    '''parse a barcode mapping file into a dictionary {barcode: sample}'''
    barcode_map = {}
    for line in barcode_lines:
        sample, barcode = line.split()
        barcode_map[barcode] = sample

    return barcode_map

def best_barcode_match(known_barcodes, barcode):
    '''
    Find the best match between a known barcode a list of known barcodes

    Parameters
    known_barcodes : sequence of iterator of sequences
        list of known barcodes
    barcode : string
        the barcode read to be matched against the known barcodes

    Returns
    min_mismatches : int
        number of mismatches in the best alignment
    best_known_barcode : string
        known barcode that aligned best
    '''
    
    # get a list of pairs (n_mismatches, known_barcode)
    n_mismatches = lambda known_barcode: mismatches(barcode, known_barcode, 1)[1]

    alignments = [(n_mismatches(known_barcode), known_barcode) for known_barcode in known_barcodes]

    # find the alignment that has the minimum number of mismatches
    min_mismatches, best_known_barcode = min(alignments, key=lambda x: x[0])

    return min_mismatches, best_known_barcode

# swo> this code is duplicated in remove_primers
def rename_fastq_ids_with_sample_names(fastq_lines, barcode_map, max_barcode_diffs):
    '''
    Rename the read IDs in a fastq file with the corresponding sample name. Get the barcode
    read right from the ID line, look it up in the barcode map, and pick the best match.



    Parameters
    fastq_lines : sequence or iterator of strings
        the lines in the fastq file to be processed
    barcode_map : dictionary
        entries are {barcode: sample_name}
    max_barcode_diffs : int
        maximum number of mismatches between a barcode read and known barcode before throwing
        out that read

    Returns nothing. Everything is printed.
    '''

    # keep track of the computations where we align the barcode read to the known barcodes
    barcode_read_to_sample = {}
    sample_counts = {}

    # get the fastq lines four at a time
    for at_line, seq_line, plus_line, quality_line in itertools.izip(*[iter(fastq_lines)] * 4):
        at_line = at_line.rstrip()
        seq_line = seq_line.rstrip()
        plus_line = plus_line.rstrip()
        quality_line = quality_line.rstrip()

        # check that the two lines with identifiers match our expectations
        assert(at_line.startswith('@'))
        assert(plus_line.startswith('+'))
        assert(len(seq_line) == len(quality_line))

        # look for the barcode from the read ID line
        # match, e.g. @any_set_of_chars#ACGT/1 -> ACGT
        m = re.match("@.*#([ACGTN]+)/\d+", at_line)

        if m is None:
            raise RuntimeError("couldn't find barcode in fastq line: %s" %(at_line))

        # if we've already aligned this barcode read, just use the same sample we found before.
        # otherwise, look through all the barcodes for the best match.
        barcode_read = m.group(1)
        if barcode_read in barcode_read_to_sample:
            sample = barcode_read_to_sample[barcode_read]
            sample_counts[sample] += 1
        else:
            # try aligning to every known barcode
            n_mismatches, best_known_barcode = best_barcode_match(barcode_map.keys(), barcode_read)

            # if the match was good, assign that barcode read to the sample that the best read
            # matches
            if n_mismatches > max_barcode_diffs:
                continue
            else:
                # get the name for this sample; record which sample we mapped this barcode
                # read to
                sample = barcode_map[best_known_barcode]
                barcode_read_to_sample[barcode_read] = sample
                sample_counts[sample] = 1

        print "@sample=%s;%d" %(sample, sample_counts[sample])
        print seq_line
        print "+"
        print quality_line


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help='input fastq file')
    parser.add_argument('barcode', help='barcode mapping file')
    parser.add_argument('-m', '--max_barcode_diffs', default=0, type=int, help='maximum number of nucleotide mismatches in the barcode (default: 0)')
    args = parser.parse_args()

    # parse the barcode mapping file
    with open(args.barcode, 'r') as f:
        barcode_map = barcode_file_to_dictionary(f)

    # get a set of reads
    with open(args.fastq, 'r') as f:
        rename_fastq_ids_with_sample_names(f, barcode_map, args.max_barcode_diffs)