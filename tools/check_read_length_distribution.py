#!/usr/bin/env python

'''
Read a fastq and get the distribution of read lengths.
'''

import argparse, sys, os
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from SmileTrain import util
from Bio import SeqIO


def fastq_count_dist(fastq):
    '''get a dictionary {length => # of reads with this length}'''
    dist = {}
    for record in SeqIO.parse(fastq, 'fastq'):
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
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help='input fastq')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output log (default stdout)')
    args = parser.parse_args()

    count_dist = fastq_count_dist(args.fastq)
    out = dist_message(count_dist)
    args.output.write(out)