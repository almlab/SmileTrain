import re, string, sys, time, itertools, os, subprocess
import usearch_python.primer

def parse_index_line(line):
    '''read an index line to sample, sequence, abundance'''
    sample, seq, abund = line.split()
    abund = int(abund)
    return [sample, seq, abund]

def parse_seq_sid(sid):
    '''seq0;counts=400 -> seq0'''
    
    m = re.match('(.*);counts=\d+', sid)
    if m is None:
        raise RuntimeError("sequence id did not parse: %s" % sid)
    else:
        return m.group(1)
    
def sid_to_sample(sid):
    '''sample=donor1;400 -> donor1'''
    
    m = re.match('sample=(.+);\d+', sid)
    if m is None:
        raise RuntimeError("fasta at line did not parse: %s" % sid)
    else:
        return m.group(1)

def counts(table, sample, otu):
    '''get counts from table structure, or 0 if sample or otu missing'''
    if sample in table:
        if otu in table[sample]:
            return table[sample][otu]
        
    return 0