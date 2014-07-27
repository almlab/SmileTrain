#!/usr/bin/env python

'''
Count the number of differences between the first sequence in a fasta file and the
following sequences. If the two sequences are not the same length, you'll get a fractional
answer.
'''

from SmileTrain import util
import argparse, difflib

def count_diffs(target, query):
    l = len(target)
    a = difflib.SequenceMatcher(a=target, b=query, autojunk=False)
    r = a.ratio()
    return l*(1.0 - r)



if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta')
    args = parser.parse_args()
    
    with open(args.fasta) as f:
        entries = list(util.fasta_entries(f))
    
    target_label, target_seq = entries[0]
    query_labels, query_seqs = zip(*entries)
    diffs = [count_diffs(target_seq, query_seq) for query_seq in query_seqs]
    max_diff = max(diffs)
    
    for query_label, diff in zip(query_labels, diffs):
        print "{}\t{}".format(query_label, diff)
        
    print "max\t{}".format(max_diff)