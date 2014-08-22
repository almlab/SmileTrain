#!/usr/bin/env python

'''
Count the number of differences between the first sequence in a fasta file and the
following sequences. If the two sequences are not the same length, you'll get a fractional
answer.
'''

from Bio import SeqIO
import argparse, difflib


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta')
    args = parser.parse_args()
    
    target = None
    max_fraction = 0.0
    for record in SeqIO.parse(args.fasta, 'fasta'):
        if target is None:
            target = str(record.seq)
            l = len(target)
        
        d = difflib.SequenceMatcher(a=target, b=str(record.seq), autojunk=False)
        fraction = l*(1.0 - d.ratio())
        print "{}\t{}".format(record.id, fraction)
        
        max_fraction = max(fraction, max_fraction)
    
    print "max\t{}".format(max_fraction)