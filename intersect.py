'''
Compare the forward and reverse fastq files. Find pairs of reads that occur in both files.
Output those reads.
'''

import sys

def fastq_entries(lines, strip=True, check_sigils=True, check_lengths=True, sequence=True):
    '''
    Reads lines from a fastq file and record them record-by-record.

    Parameters
    lines : sequence or iterator of strings
        lines from the fastq file
    strip : bool (default: true)
        rstrip all the lines?
    check_sigils : bool (default: true)
        assert that the @ and + lines begin with those symbols?
    check_lengths : bool (default: true)
        assert that the sequence and quality lines are the same length?
    sequence : bool (default: true)
        yield a sequence of four strings; otherwise, a single newline-separated string

    iterator
    '''

    for at_line, seq_line, plus_line, quality_line in itertools.izip(*[iter(lines)] * 4):
        if strip:
            # chomp all newlines
            at_line = at_line.rstrip()
            seq_line = seq_line.rstrip()
            plus_line = plus_line.rstrip()
            quality_line = quality_line.rstrip()

        if check_sigils:
            # check that the two lines with identifiers match our expectations
            assert(at_line.startswith('@'))
            assert(plus_line.startswith('+'))

        if check_lengths:
            # check that the sequence and quality lines have the same number of nucleotides
            assert(len(seq_line) == len(quality_line))

        # return a list or a single string
        if sequence:
            yield [at_line, seq_line, plus_line, quality_line]
        else:
            yield "\n".join([at_line, seq_line, plus_line, quality_line])

def fastq_ids(lines):
    '''Extract the read IDs (with no @) from a fastq file'''
    return [at_line[1:] for at_line, seq_line, plus_line, quality_line in fastq_entries(lines)]

def common_ids(fastq1, fastq2):
    '''
    Get a list of IDs corresponding to reads found in the two fastq files.

    Parameters
    fastq1, fastq2 : sequence or iterator of strings
        lines from the fastq files

    Returns
    common_ids : set
        ids (not including the @)
    '''

    ids1 = fastq_ids(fastq1)
    ids2 = fastq_ids(fastq2)

    common_ids = set([id1 for id1 in ids1 if id1 in ids2])
    return common_ids

def fastq_entries_with_matching_ids(fastq, rids):
    '''
    Yield a series of fastq entries drawn from the input whose IDs match those in the list.
    Each entry is a single string.
    '''

    for at_line, seq_line, plus_line, quality_line in fastq_entries(fastq):
        rid = at_line[1:]

        if rid in rids:
            yield "\n".join([at_line, seq_line, plus_line, quality_line]) 


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('forward_in', help='input forward reads fastq file')
    parser.add_argument('reverse_in', help='input reverse reads fastq file')
    parser.add_argument('forward_out', help='output forward reads fastq file')
    parser.add_argument('reverse_out', help='output reverse reads fastq file')
    args = parser.parse_args()

    # read in the inputs and look for common ids
    with open(args.forward_in, 'r') as f:
        with open(args.reverse_in, 'r') as r:
            rids = common_ids(f, r)

    # write the forward entries with reads with ids in the reverse entries
    with open(args.forward_in, 'r') as i:
        with open(args.forward_out, 'w') as o:
            for entry in fastq_entries_with_matching_ids(i, rids):
                o.write(entry)

    # ditto for the revere reads
    with open(args.reverse_in, 'r') as i:
        with open(args.reverse_out, 'w') as o:
            for entry in fastq_entries_with_matching_ids(i, rids):
                o.write(entry)