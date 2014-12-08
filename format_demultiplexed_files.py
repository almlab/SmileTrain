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
    
and then the resulting fastq files are concatenated and placed in the output file

If the barcode map has three fields, corresponding to forward, reverse, and sample

    filename1.F.fastq filename1.R.fastq sample_a
    filename2.F.fastq filename2.R.fastq sample_b

the forward and reverse files will be concatenated separately and output appropriately
'''

import usearch_python.primer, util
import sys, argparse, string, itertools, re
from Bio import SeqIO

def barcode_file_to_dictionary(barcode_lines):
    '''parse a filename mapping file into a dictionary {filename: sample}'''
    barcode_map = {}
    barcode_map_forward = {}
    barcode_map_reverse = {}

    for i, line in enumerate(barcode_lines):
        fields = line.split()

        if len(fields) == 2:
            filename, sample = fields
            barcode_map[filename] = sample

        if len(fields) == 3:
            forwardFile, reverseFile, sample = fields
            barcode_map_forward[forwardFile] = sample
            barcode_map_reverse[reverseFile] = sample

        if len(fields) not in (2,3):
            raise RuntimeError('All lines in mapping file must have 2 or three fields, found %d in line: %s' %(len(fields, line)))

    if barcode_map_forward and barcode_map_reverse:
        return barcode_map_forward, barcode_map_reverse
    else:
        return barcode_map

def renamed_fastq_records(fastq, barcode_map, reverse=False):
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
    count = 0
    for record in SeqIO.parse(fastq, 'fastq'):

        if not reverse:
            desc = '%d#%s/1' %(count, barcode_map[fastq])
        else:
            desc = '%d#%s/2' %(count, barcode_map[fastq])
        record.id = desc
        record.description = ''
        count += 1
        yield record

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reformat and concatenate pre-demultiplexed files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('barcode', help='barcode mapping file')
    parser.add_argument('--output', '-o', default=None, help='output concatenated fastq')
    args = parser.parse_args()

    # parse the barcode mapping file
    with open(args.barcode, 'r') as f:
        barcode_map = barcode_file_to_dictionary(f)

    if type(barcode_map) is dict:
        # get a set of reads
        with open(args.output, 'w') as outputFile:
            for filename in barcode_map.keys():
                for record in renamed_fastq_records(filename, barcode_map):
                    SeqIO.write(record, outputFile, 'fastq')
    
    if type(barcode_map) is tuple:
        barcode_map_forward, barcode_map_reverse = barcode_map
        with open('f.'+args.output, 'w') as outputFile:
            for filename in barcode_map_forward.keys():
                for record in renamed_fastq_records(filename, barcode_map_forward):
                    SeqIO.write(record, outputFile, 'fastq')
        with open('r.'+args.output, 'w') as outputFile:
            for filename in barcode_map_reverse.keys():
                for record in renamed_fastq_records(filename, barcode_map_reverse, reverse=True):
                    SeqIO.write(record, outputFile, 'fastq')
