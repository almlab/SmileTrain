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

import sys

def fasta_entries(lines):
    '''
    Yield [id, sequence] pairs from a fasta file. Sequence allowed to run over multiple
    lines.

    Parameters
    lines : sequence or iterator of strings
        lines from the fasta file

    Yields [id (without the >), sequence] pairs
    '''

    entry = ''
    for line in lines:
        line = line.rstrip()
        if line.startswith('>'):
            if entry != '':
                yield [sid, sequence]
            else:
                sid = line[1:]
                sequence = ''
        else:
            sequence += line

    yield [sid, sequence]

def sequence_abundances(entries):
    '''fasta entries to dictionary {sequence: # reads}'''
    abundance_map = {}
    for sid, sequence in entries:
        if sequence not in abundance_map:
            abundance_map[sequence] = 1
        else:
            abundance_map[sequence] += 1

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


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input fasta file')
    parser.add_argument('output', help='output fasta file')
    parser.add_argument('-m', '--minimum_counts', type=int, default=2, help='minimum times a sequence is included, otherwise it gets thrown out (default: 2)')
    args = parser.parse_args()

    with open(args.input, 'r') as f:
        abundance_map = sequence_abundances(fasta_entries(f))

    sorted_seqs = sorted_abundant_keys(abundance_map, args.minimum_counts)

    with open(args.output, 'w') as f:
        for otu_i, seq in enumerate(sorted_seqs):
            f.write('>otu%s;size=%d\n' % (otu_i, abundance_map[seq]))
            f.write('%s\n' %(seq))