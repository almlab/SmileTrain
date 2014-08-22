#!/usr/bin/env python

'''
some index fastq's have a weird number of quality line characters. some have an extra
character; others seem to have a single character.

this script truncates quality lines longer than the sequence line and pads quality
lines that are shorter than the sequence line.

author : scott w olesen <swo@mit.edu>
'''

import argparse, sys, os, itertools
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser(description='correct quality line length')
    parser.add_argument('fastq', help='input barcode fastq')
    parser.add_argument('-z', '--fill_char', default='F', help='fill character (default: F)')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'), help='output fastq (default: stdout)')
    args = parser.parse_args()
    
    with open(args.fastq) as f:
        for four_lines in itertools.izip(*[iter(f)]*4):
            at_line, seq_line, plus_line, quality_line = [l.rstrip() for l in four_lines]
            ls = len(seq_line)
            lq = len(quality_line)
            if lq < ls:
                quality_line = quality_line.ljust(len(seq_line), args.fill_char)
            elif lq > ls:
                quality_line = quality_line[0: ls]
        
            args.output.write("\n".join([at_line, seq_line, plus_line, quality_line]) + "\n")