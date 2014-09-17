#!/usr/bin/env python

'''
Read a fastq and output the distribution of quality scores.
'''

import argparse, sys, os
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from SmileTrain import util
from Bio import SeqIO

import numpy as np

def boxplot_stats(a):
    out = {}
    out['Q1'], out['Q2'], out['Q3'] = np.percentile(a, [25, 50, 75])
    iqr = out['Q3'] - out['Q1']

    candidate_whiskers = np.logical_and(a <= out['Q1'], a >= out['Q1'] - 1.5*iqr)
    if np.any(candidate_whiskers):
        out['lower whisker'] = np.min(a[candidate_whiskers])
    else:
        out['lower whisker'] = out['Q1']

    candidate_whiskers = np.logical_and(a >= out['Q3'], a <= out['Q3'] + 1.5*iqr)
    if np.any(candidate_whiskers):
        out['upper whisker'] = np.min(a[candidate_whiskers])
    else:
        out['upper whisker'] = out['Q3']

    out['outliers'] = a[np.logical_or(a < out['lower whisker'], a > out['upper whisker'])]

    return out

def fastq_quality_dist(fastq):
    record_qualities = [record.letter_annotations['phred_quality'] for record in SeqIO.parse(fastq, 'fastq')]

    # get the maximum read length
    max_len = max([len(x) for x in record_qualities])

    # pad the qualities with nan's
    for q in record_qualities:
        q += [np.nan] * (max_len - len(q))

    position_qualities = np.vstack(tuple(record_qualities)).T
    position_stats = [boxplot_stats(position) for position in position_qualities]

    return position_stats

def boxplot_line(stats, plot_min=33, plot_max=73):
    out = "|"
    for i in range(plot_min, plot_max+1):
        if i < stats['lower whisker']:
            out += " "
        elif i >= stats['lower whisker'] and i < stats['Q1']:
            out += "-"
        elif i >= stats['Q1'] and i < stats['Q2']:
            out += "="
        elif i == stats['Q2']:
            out += "0"
        elif i > stats['Q2'] and i <= stats['Q3']:
            out += "="
        elif i > stats['Q3'] and i <= stats['upper whisker']:
            out += "-"
        elif i > stats['upper whisker']:
            out += " "

    out += "|"
    return out

def boxplot(position_stats):
    lines = [quality_header()]
    lines += ["%s %d" %(boxplot_line(stat), i+1) for i, stat in enumerate(position_stats)]
    return "\n".join(lines) + "\n"

def quality_header():
    return "|" + "".join([chr(x+1) for x in range(33, 73+1)]) + "|"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="get a fastq quality report")
    parser.add_argument('fastq', help='input fastq')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output log (default stdout)')
    args = parser.parse_args()

    stats = fastq_quality_dist(args.fastq)
    print boxplot(stats)