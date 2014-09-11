#!/usr/bin/env python

'''
Compare the forward and reverse fastq files. Make sure that there are pairs of forward and
reverse reads that are appropriately labeled.
'''

import sys, argparse, itertools
from Bio import SeqIO
import util

def check_fastq_id(fid, direction):
    if not fid.endswith("/%s" %(direction)):
        raise RuntimeError("fastq id not correct direction: '%s' should end in '%s'" %(fid, direction))

def check_for_paired_ids(forward_fastq, reverse_fastq):
    for forward_record, reverse_record in itertools.izip(SeqIO.parse(forward_fastq, 'fastq'), SeqIO.parse(reverse_fastq, 'fastq')):
        check_fastq_id(forward_record.id, '1')
        check_fastq_id(reverse_record.id, '2')


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('forward', help='forward reads fastq')
    parser.add_argument('reverse', help='reverse reads fastq')
    args = parser.parse_args()

    check_for_paired_ids(args.forward, args.reverse)