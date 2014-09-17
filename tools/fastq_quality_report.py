#!/usr/bin/env python

'''
Read a fastq and output the distribution of quality scores.
'''

import argparse, sys, os
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from SmileTrain import util
from Bio import SeqIO

def quality_dictionary():
    '''container for counting quality scores {33 => 0, ..., 73 => 0}'''
    return {i: 0 for i in range(33, 73+1)}

def fastq_quality_dist(fastq):
    positions = []

    for record in SeqIO.parse(fastq, 'fastq'):
        qualities = record.letter_annotations['phred_quality']

        # add on empty quality containers until we match the read length
        while len(qualities) > len(positions):
            positions.append(quality_dictionary())

        for position, quality in zip(positions, qualities):
            position[quality] += 1

    return positions


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="get a fastq quality report")
    parser.add_argument('fastq', help='input fastq')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output log (default stdout)')
    args = parser.parse_args()

    positions = fastq_quality_dist(args.fastq)
    print positions