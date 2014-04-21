'''
Compare the forward and reverse fastq files. Find pairs of reads that occur in both files.
Output those reads.
'''

import sys
import util

def fastq_ids(lines):
    '''Extract the read IDs (with no @) from a fastq file'''
    return [at_line[1:] for at_line, seq_line, quality_line in util.fastq_iterator(lines)]

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

    ids1 = set(fastq_ids(fastq1))
    ids2 = set(fastq_ids(fastq2))

    return ids1.intersection(ids2)

def fastq_entries_with_matching_ids(fastq, rids):
    '''yield a series of fastq entry strings drawn from the input whose IDs match those in the list'''

    for at_line, plus_line, quality_line in util.fastq_iterator(fastq):
        rid = at_line[1:]

        if rid in rids:
            yield util.fastq_entry_list_to_string([at_line, seq_line, quality_line]) 


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
