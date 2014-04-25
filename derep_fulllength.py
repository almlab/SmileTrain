'''
Dereplicate a fasta file. Read through all the entries, finding reads with identical
sequences. Write out a new fasta file with IDs that keep track of the abundance of
each sequence.

input like
    >foo
    CATCAT
    >bar
    TATA
    >baz
    AGGAGG
    >lorem
    TATA
    >ipsum
    CATCAT
    >qualcum
    CATCAT

with a minimum number of reads
    2

becomes output like
    >otu1;size=3
    CATCAT
    >otu2;size=2
    TATA

The output file is sorted in order of decreasing abundance. Sequences with a minimum
number of counts are dropped.
'''

import sys, argparse, re
import util

def sid_to_sample(sid):
    '''>sample=donor1;400 -> donor1'''
    
    m = re.match('sample=(.+);\d+', sid)
    if m is None:
        raise RuntimeError("fasta at line did not parse: %s" % sid)
    else:
        return m.group(1)

def sequence_abundances(entries, by_sample=False):
    '''
    fasta entries to dictionary {sequence: # reads}
    
    
    '''
    
    abundance_map = {}
    for sid, sequence in entries:
        if by_sample:
            key = (sid_to_sample(sid), sequence)
        else:
            key = sequence
        
        if key not in abundance_map:
            abundance_map[key] = 1
        else:
            abundance_map[key] += 1

    return abundance_map

def sorted_abundant_keys(abundance_map, minimum_counts):
    '''
    Given a dictionary {key: number}, get a list of keys. The keys are sorted so that
    d[key] are numbers in decreasing order and d[key] is at least minimum counts.

    Parameters
    abundance_map : dictionary
        a dictionary {key: number}
    minimum_counts : int
        minimum number to avoid being dropped

    Returns
    keys : list
        keys in decreasing order with numbers above the minimum
    '''

    sorted_keys = sorted(abundance_map.keys(), key=lambda k: abundance_map[k], reverse=True)
    filtered_sorted_keys = [k for k in sorted_keys if abundance_map[k] >= minimum_counts]

    return filtered_sorted_keys

def sorted_abundant_entries(entries, min_counts, by_sample=False):
    '''
    Given fasta entries, find the most abundant ones, and return them in order.
    
    entries : list or iterator of [name, sequence] pairs
        fasta entries
    min_counts : int
        minimum number of times a sequence must appear to be counted
    by_sample : bool (default false)
        only dereplicate sequences in the same sample
        
    yields : lists
        [new name, sequence] pairs
    '''
    
    abundance_map = sequence_abundances(entries, by_sample=by_sample)
    sorted_keys = sorted_abundant_keys(abundance_map, min_counts)
    
    if by_sample:
        seq_i = {}
        for key in sorted_keys:
            sample, seq = key
            if sample not in seq_i:
                seq_i[sample] = 1
            else:
                seq_i[sample] += 1
                
            yield ">sample=%s;seq=%d;counts=%d\n%s" %(sample, seq_i[sample], abundance_map[key], seq)
    else:
        for otu_i, seq in enumerate(sorted_keys):
            yield '>otu%s;counts=%d\n%s' %(otu_i, abundance_map[seq], seq)


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input fasta file')
    parser.add_argument('output', help='output fasta file')
    parser.add_argument('-m', '--minimum_counts', type=int, default=2, help='minimum times a sequence is included, otherwise it gets thrown out (default: 2)')
    parser.add_argument('--by_sample', action='store_true', help='only dereplicate the same sequence when found in the same sample')
    args = parser.parse_args()

    with open(args.input, 'r') as f:
        with open(args.output, 'w') as o:
            for entry in sorted_abundant_entries(util.fasta_entries(f), args.minimum_counts, by_sample=args.by_sample):
                o.write("%s\n" %(entry))