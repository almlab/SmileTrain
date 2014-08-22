#!/usr/bin/env python

'''
Compare the forward and reverse fastq files. Find pairs of reads that occur in both files.
Output those reads.

All @ lines are being parsed in the same way. Depending on the file types that we have
being fed into the pipeline, that naming convention would need to change.
'''

import sys, argparse, re
from Bio import SeqIO
import util

def fastq_id_to_read_id(fid):
    '''trim off the /1 or /2'''
    m = re.search("/[12]$", fid)
    if not m:
        raise RuntimeError("could not parse fastq id: %s" %(fid))
    
    return fid.rstrip("/12")

def fastq_ids(fastq):
    '''extract the read IDs from a fastq file'''
    return [fastq_id_to_read_id(record.id) for record in SeqIO.parse(fastq, 'fastq')]

def common_ids(fastq1, fastq2):
    '''
    Get a list of IDs corresponding to reads found in the two fastq files.

    Parameters
    fastq1, fastq2 : sequence or iterator of strings
        lines from the fastq files

    Returns
    common_ids : set
        ids (not including the @)
    '''

    ids1 = set(fastq_ids(fastq1))
    ids2 = set(fastq_ids(fastq2))

    return ids1.intersection(ids2)

def fastq_records_with_matching_ids(fastq, rids):
    '''yield a series of fastq entry strings drawn from the input whose IDs match those in the list'''
    
    for record in SeqIO.parse(fastq, 'fastq'):
        if fastq_id_to_read_id(record.id) in rids:
            yield record


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('forward_in', help='input forward reads fastq file')
    parser.add_argument('reverse_in', help='input reverse reads fastq file')
    parser.add_argument('forward_out', help='output forward reads fastq file')
    parser.add_argument('reverse_out', help='output reverse reads fastq file')
    args = parser.parse_args()

    # read in the inputs and look for common ids
    rids = common_ids(args.forward_in, args.reverse_in)
            
    # make sure that we are actually looking for something
    if len(rids) == 0:
        raise RuntimeError("no common IDs found in %s and %s" % (args.forward_in, args.reverse_in))

    # write the forward entries with reads with ids in the reverse entries
    with open(args.forward_out, 'w') as f:
        for record in fastq_records_with_matching_ids(args.forward_in, rids):
            SeqIO.write(record, f, 'fastq')
        
    with open(args.reverse_out, 'w') as f:
        for record in fastq_records_with_matching_ids(args.reverse_in, rids):
            SeqIO.write(record, f, 'fastq')