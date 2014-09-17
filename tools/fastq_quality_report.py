#!/usr/bin/env python

'''
Read a fastq and output the distribution of quality scores.
'''

import argparse, sys, os, textwrap
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from SmileTrain import util
from Bio import SeqIO

import numpy as np

# range of possible read quality scores
min_q = 1
max_q = 41

def boxplot_stats(a):
    '''given an array, return a dictionary with some stats useful for box plotting'''
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

def empty_position():
    '''empty dictionary for holding quality score counts'''
    return {x: 0 for x in range(min_q, max_q+1)}

def quality_dictionary_to_array(d):
    '''{1 => 3, 2 => 1} -> [1, 1, 1, 2]'''
    o = []
    for quality, count in d.items():
        o += [quality]*count

    return np.array(o)

def fastq_quality_dist(fastq):
    '''
    Read a fastq to create a list of dictionaries. The first dictionaries is {quality => counts}
    for the first position in the reads, etc.
    '''

    position_qualities = []

    for record in SeqIO.parse(fastq, 'fastq'):
        quality = record.letter_annotations['phred_quality']

        while len(quality) > len(position_qualities):
            position_qualities.append(empty_position())

        for q, position in zip(quality, position_qualities):
            if q > max_q or q < min_q:
                raise RuntimeError("Quality %d outside of expected range (%d, %d)" %(q, min_q, max_q))

            position[q] += 1

    return position_qualities

def fastq_quality_stats(fastq):
    '''
    Read a fastq to create a list of stats dictionaries. Every item in the list is for one
    position amongst the reads (e.g., first item for first bp across all reads).
    '''

    position_qualities = fastq_quality_dist(fastq)
    position_stats = [boxplot_stats(quality_dictionary_to_array(pq)) for pq in position_qualities]

    return position_stats


def boxplot_line(stats, plot_min=min_q, plot_max=max_q):
    '''given statistics, write out a horizontal line in the boxplot'''
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
    '''given stats for each position, make out a horizontal boxplot'''
    lines = [quality_header()]
    lines += ["%s %d" %(boxplot_line(stat), i+1) for i, stat in enumerate(position_stats)]
    return "\n".join(lines) + "\n"

def quality_header():
    '''boxplot header that shows quality score characters'''
    return "|" + "".join([chr(x+1) for x in range(33, 73+1)]) + "|"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            Get a fastq quality report. Every line is one position in the reads.

            0 median
            = within IQR
            - within 1.5 IQR of Q1 or Q3
            '''))
    parser.add_argument('fastq', help='input fastq')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output log (default stdout)')
    args = parser.parse_args()

    stats = fastq_quality_stats(args.fastq)
    print boxplot(stats)