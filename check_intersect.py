#!/usr/bin/env python

'''
Compare the forward and reverse fastq files. Make sure that there are pairs of forward and
reverse reads that are appropriately labeled.
'''

import sys, argparse, itertools, subprocess
from Bio import SeqIO
import util

def n_lines(fn):
    '''count number of lines in fn using wc'''
    return int(subprocess.check_output(['wc', '-l', fn]).split()[0])

def check_for_matched_line_numbers(f1, f2):
    '''assert that two files have the same number of lines'''
    if n_lines(f1) != n_lines(f2):
        raise RuntimeError("fastqs have unequal numbers of lines")

def check_matching_fastq_ids(forward_id, reverse_id):
    '''check that reads are whatever/1 and whatever/2'''
    if not forward_id.endswith("/1"):
        raise RuntimeError("fastq id not correct direction: '%s' should end in 1" % forward_id)

    if not reverse_id.endswith("/2"):
        raise RuntimeError("fastq id not correct direction: '%s' should end in 2" % reverse_id)

    if forward_id.rstrip("/12") != reverse_id.rstrip("/12"):
        raise RuntimeError("fastq ids did not match: '%s' and '%s'" % (forward_id, reverse_id))

def check_for_paired_ids(forward_fastq, reverse_fastq):
    '''check each pair of fastq records for paired read ids'''
    for forward_record, reverse_record in itertools.izip(SeqIO.parse(forward_fastq, 'fastq'), SeqIO.parse(reverse_fastq, 'fastq')):
        check_matching_fastq_ids(forward_record.id, reverse_record.id)


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('forward', help='forward reads fastq')
    parser.add_argument('reverse', help='reverse reads fastq')
    args = parser.parse_args()

    check_for_matched_line_numbers(args.forward, args.reverse)
    check_for_paired_ids(args.forward, args.reverse)