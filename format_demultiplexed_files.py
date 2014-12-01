#!/usr/bin/env python

'''
Reformat header of demultiplexed files and concatenate them

Given a tab-separated filename mapping file like
    filename1.fastq   sample_a
    filename2.fastq   sample_b

the headers of sequences in filename1.fastq:
    @OURSEQ:lolapalooza1234/1
    AACCGGTT
    +
    abcdefgh

becomes output like
    @OURSEQ:lolapalooza1234#sample_a/1
    AACCGGTT
    +whatever
    abcdefgh
    
and then the resulting fastq files are concatenated and placed int he output file
'''

import usearch_python.primer, util
import sys, argparse, string, itertools, re
from Bio import SeqIO

def barcode_file_to_dictionary(barcode_lines):
    '''parse a filename mapping file into a dictionary {filename: sample}'''
    barcode_map = {}
    for i, line in enumerate(barcode_lines):
        fields = line.split()

        if len(fields) != 2:
            raise RuntimeError("every line in filename mapping file should have two fields; found %d in line %d" %(len(fields), i))

        filename, sample = fields
        barcode_map[filename] = sample

    return barcode_map

def renamed_fastq_records(fastq, barcode_map):
    '''
    Rename the read IDs in a fastq file with the corresponding sample name.

    Parameters
    fastq : filename or filehandle
        input
    barcode_map : dictionary
        entries are {filename: sample}

    yields : SeqRecord
        fastq records
    '''

    whitespace = re.compile(r'\s+')

    # keep track of the computations where we align the barcode read to the known barcodes

    for record in SeqIO.parse(fastq, 'fastq'):

        desc = whitespace.sub(':',record.description)
        desc += '#%s/1' %(barcode_map[fastq])
        record.id = desc
        record.description = ''
        yield record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reformat and concatenate pre-demultiplexed files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('barcode', help='barcode mapping file')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output concatenated fastq')
    args = parser.parse_args()

    # parse the barcode mapping file
    with open(args.barcode, 'r') as f:
        barcode_map = barcode_file_to_dictionary(f)

    # get a set of reads
    for filename in barcode_map.keys():
        for record in renamed_fastq_records(filename, barcode_map):
            SeqIO.write(record, args.output, 'fastq')