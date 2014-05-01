#!/usr/bin/env python

'''
Create an OTU table by combining information from
    * the .uc file that contains the mapping sequence ID => OTU
    * the index file with entries (sample, sequence ID, counts)
'''

import re, sys, argparse

def parse_uc_line(line):
    '''uc line -> (sequence ID, OTU)'''
    fields = line.split()
    
    label = fields[8]
    
    m = re.match('(.*);counts=\d+', label)
    if m is None:
        raise RuntimeError("uc label did not parse: %s" % label)
    else:
        sid = m.group(1)
    
    hit = fields[0]
    if hit == 'H':
        otu = fields[9]
    elif hit == 'N':
        otu = 'no_match'
    else:
        raise RuntimeError('unknown code %s found in .uc file' % hit)
    
    return (sid, otu)

def parse_uc_lines(lines):
    '''uc lines -> dictionary {sequence ID => OTU}'''
    
    sid_otu = {}
    for line in lines:
        sid, otu = parse_uc_line(line)
        
        if sid in sid_otu:
            raise RuntimeError("sequence ID %s repeated in .uc lines" % sid)

        sid_otu[sid] = otu
    
    return sid_otu

def parse_sample_lines(lines):
    '''read in the first column of each line'''
    return [line.split()[0] for line in lines]

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
        sample, seq, abund = parse_index_line(line)
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
        yield "\t".join([otu] + [str(counts(table, sample, otu)) for sample in samples])


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
