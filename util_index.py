import re, string, sys, time, itertools, os, subprocess
import usearch_python.primer

def parse_index_line(line):
    '''read an index line to sample, sequence, abundance'''
    sample, seq, abund = line.split()
    abund = int(abund)
    return [sample, seq, abund]

def counts(table, sample, otu):
    '''get counts from table structure, or 0 if sample or otu missing'''
    if sample in table:
        if otu in table[sample]:
            return table[sample][otu]
        
    return 0