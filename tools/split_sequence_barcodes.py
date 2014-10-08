#!/usr/bin/env python

'''
Search for each of the barcodes in the sequences. Send sequences mapping different barcodes
into different fastq files.
'''

from SmileTrain.usearch_python import primer as usearch_primer
from SmileTrain import util, map_barcodes
import sys, argparse, string, itertools, re
from Bio import SeqIO, SeqRecord

def best_plate_barcode_match(known_barcodes, seq, window=15):
    '''
    Find the best match between a known barcode a list of known barcodes

    Parameters
    known_barcodes : sequence of iterator of sequences
        list of known barcodes
    seq : string
        the sequence read to be matched against the known barcodes
    window : int
        how far into the seq should you look for the barcode

    Returns
    best_barcode

    best_pos

    best_diffs

    '''

    alignments = [(bc, pos, usearch_primer.MatchPrefix(seq, bc[pos:])) for pos in range(window) for bc in known_barcodes]

    best_barcode, best_pos, best_diffs = min(alignments, key=lambda x: x[2])

    return best_barcode, best_pos, best_diffs

def split_records(fastq, barcode_map, max_barcode_diffs):
    '''
    Split the records into new files based on their plate barcodes.

    Parameters
    fastq : filename or filehandle
        input
    barcode_map : dictionary
        entries are {barcode: destination filehandle or filename}
    max_barcode_diffs : int
        maximum number of mismatches between a barcode read and known barcode before throwing
        out that read

    yields : SeqRecord
        fastq records
    '''

    for record in SeqIO.parse(fastq, 'fastq'):
        # look for the barcode from the read ID line
        seq = str(record.seq)
        known_barcodes = barcode_map.keys()
        best_barcode, best_position, best_diffs = best_plate_barcode_match(known_barcodes, seq)
        
        # if the match was good, assign that barcode read to the sample that the best read
        # matches
        if best_diffs <= max_barcode_diffs:
            quality = record.letter_annotations['phred_quality']

            barcode_end_index = best_position + len(best_barcode)
            trim_seq = seq[barcode_end_index:]
            trim_quality = quality[barcode_end_index:]
            new_record = SeqRecord.SeqRecord(seq=str(trim_seq), id=record.id, letter_annotations={'phred_quality': trim_quality}, description="")

            destination = barcode_map[best_barcode]

            SeqIO.write(new_record, destination, 'fastq')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Demultiplex fastq entries by plate barcode', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fastq', help='input fastq file')
    parser.add_argument('barcode', help='barcode mapping file')
    parser.add_argument('-m', '--max_barcode_diffs', default=0, type=int, help='maximum number of nucleotide mismatches in the barcode')
    args = parser.parse_args()

    # parse the barcode mapping file
    with open(args.barcode, 'r') as f:
        name_barcode_map = map_barcodes.barcode_file_to_dictionary(f)

    barcode_map = {bc: open("%s.%s" % (args.fastq, name_barcode_map[bc]), 'w') for bc in name_barcode_map}
    print barcode_map

    split_records(args.fastq, barcode_map, args.max_barcode_diffs)