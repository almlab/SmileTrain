#!/usr/bin/env python

'''
Read a fastq and get the distribution of read lengths.
'''

import argparse, sys, os
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from SmileTrain import util
from Bio import SeqIO


def count_dist(fn, ftype):
    '''
    Distribution of counts in fasta or fastq

    Parameters
    fn : filename or filehandle
        input file
    ftype : 'fastq' or 'fasta'
        input file type

    returns: dictionary
        {length => # of reads with this length}

    '''

    dist = {}
    for record in SeqIO.parse(fn, ftype):
        l = len(record.seq)
        if l in dist:
            dist[l] += 1
        else:
            dist[l] = 1

    return dist

def dist_message(dist):
    '''a nice output message based on the read distribution'''
    total = sum(dist.values())
    dist = sorted(list(dist.items()))

    lines = ["length\tcounts"]
    lines += ["%d\t%d" %(length, counts) for length, counts in dist]
    out = "\n".join(lines) + "\n"
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read fastq or fasta and get distribution of sequence lengths', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    type_group = parser.add_mutually_exclusive_group(required=True)
    type_group.add_argument('-a', '--fasta', help='input fasta')
    type_group.add_argument('-q', '--fastq', help='input fastq')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output log')
    args = parser.parse_args()

    if args.fasta:
        count_dist = count_dist(args.fasta, 'fasta')
    elif args.fastq:
        count_dist = count_dist(args.fastq, 'fastq')
        
    out = dist_message(count_dist)
    args.output.write(out)