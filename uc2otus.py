#!/usr/bin/env python

'''
Create an OTU table by combining information from
    * the .uc file that contains the mapping sequence ID => OTU
    * the index file with entries (sample, sequence ID, counts)
'''

import re, sys, argparse
import util_index

def parse_uc_line(line):
    '''
    uc line -> (H or N, sid, otu)
    '''
    fields = line.split()
    
    hit = fields[0]
    label = fields[8]
    otu = fields[9]
    
    return (hit, label, otu)

def parse_uc_lines(lines, miss_name='no_match'):
    '''uc lines -> dictionary {sequence ID => OTU}'''
    
    sid_otu = {}
    for line in lines:
        # rename the otu from "*" if there was no hit
        hit, label, otu = parse_uc_line(line)
        if hit == 'N':
            otu = miss_name
        elif hit != 'H':
            raise RuntimeError('unknown code %s found in .uc file' % hit)
    
        # parse the sid "seq123;counts=456" to "seq123"    
        m = re.match('(.*);(counts|size)=\d+', label)
        if m is None:
            raise RuntimeError("uc label did not parse: %s" % label)
        else:
            sid = m.group(1)
        
        if sid in sid_otu:
            raise RuntimeError("sequence ID %s repeated in .uc lines" % sid)

        sid_otu[sid] = otu
    
    return sid_otu

def parse_sample_lines(lines):
    '''read in the first column of each line'''
    return [line.split()[0] for line in lines]

def sparse_count_table(seq_otu, index_lines):
    '''
    Create a sparse count table using sequence-OTU mapping and sample-sequence-abundance index.
    
    seq_otu : dictionary
        {'sequence id' => 'otu name'}
    index_lines : list or iterator of strings
        lines from the index file (tab-separated sample, sequence ID, abundance)
        
    returns : dictionary of dictionaries
        {sample => {otu => abundance, ...}, ...}
    '''
    
    table = {}
    for line in index_lines:
        sample, seq, abund = util_index.parse_index_line(line)
        otu = seq_otu[seq]
        
        if sample not in table:
            table[sample] = {otu: abund}
        else:
            if otu not in table[sample]:
                table[sample][otu] = abund
            else:
                table[sample][otu] += abund
                
    return table

def otu_table(table, otus=None, samples=None):
    '''
    Output lines of an OTU table.
    
    table : dictionary of dictionaries
        {sample => {otu => abundance, ...}, ...}
    otus : list or iterator of strings (default None)
        otu ids in row order; or just make a new sorted list
    samples : list or iterator of strings (default None)
        sample names in column order; or just make a new sorted list
        
    yields : strings
        lines in the otu table
    '''
    
    # make our own otu list if necessary
    if samples is None:
        samples = sorted(table.keys())
        
        # if 'no_match' is in there, put it first
        if 'no_match' in samples:
            samples.remove('no_match')
            samples.insert(0, 'no_match')
        
    # and our own samples list
    if otus is None:
        # concatenate the lists of keys
        otus = []
        for sample in table:
            otus += table[sample].keys()
            
        # remove duplicates and sort
        otus = sorted(list(set(otus)))
    
    # first, output the header/sample line
    yield "\t".join(['OTU_ID'] + samples)
    
    # loop over rows
    for otu in otus: 
        yield "\t".join([otu] + [str(util_index.counts(table, sample, otu)) for sample in samples])


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('uc', help='input uc file')
    parser.add_argument('index', help='input index file')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    parser.add_argument('--samples', default=None, help='samples in order in the first field (e.g., a barcode file)')
    parser.add_argument('--otus', default=None, help='OTUs in order in the first field')
    args = parser.parse_args()
    
    with open(args.uc) as f:
        sid_otu = parse_uc_lines(f)
        
    with open(args.index) as f:
        table = sparse_count_table(sid_otu, f)
        
    if args.samples is None:
        samples = None
    else:
        with open(args.samples) as f:
            samples = parse_sample_lines(f)
    
    if args.otus is None:
        otus = None
    else:
        with open(args.otus) as f:
            otus = parse_sample_lines(f)
            
    for line in otu_table(table, otus, samples):
        args.output.write("%s\n" % line)
