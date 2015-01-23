#!/usr/bin/env python

'''
Reformat header of demultiplexed files and concatenate them

Given a tab-separated filename mapping file like
    filename1.fastq   sample_a
    filename2.fastq   sample_b

the headers of sequences in filename1.fastq:
    @OURSEQ:lolapalooza1234/1
    AACCGGTT
    +
    abcdefgh

becomes output like
    @OURSEQ:lolapalooza1234#sample_a/1
    AACCGGTT
    +whatever
    abcdefgh

and then the resulting fastq files are concatenated and placed in the output file

If the sample map has three fields, corresponding to forward, reverse, and sample

    filename1.F.fastq filename1.R.fastq sample_a
    filename2.F.fastq filename2.R.fastq sample_b

the forward and reverse files will be concatenated separately and output appropriately
'''

import argparse
from itertools import izip_longest
from SmileTrain import util
from Bio import SeqIO


def sample_file_to_dictionary(sample_file):
    '''parse a filename mapping file into a dictionary {filename: sample}'''
    sample_map = {}
    lines = sample_file.readlines()
    num_fields = len(lines[0].split())

    if num_fields == 2:
        for i, line in enumerate(lines):
            fields = line.split()
            if not len(fields) is 2:
                raise RuntimeError('All lines in mapping file must have same number of fields. First line has 2, found %d in line: %s' % (len(fields), line))
            filename, sample = fields
            sample_map[filename] = sample
        return sample_map

    elif num_fields == 3:
        for i, line in enumerate(lines):
            fields = line.split()
            if not len(fields) is 3:
                raise RuntimeError('All lines in mapping file must have same number of fields. First line has 3, found %d in line: %s' % (len(fields), line))
            forwardFile, reverseFile, sample = fields
            sample_map[(forwardFile, reverseFile)] = sample
        return sample_map

    elif num_fields not in (2, 3):
        raise RuntimeError('All lines in mapping file must have 2 or three fields, found %d on first line' % (len(fields)))


def renamed_fastq_records(key, sample_map, reverse=False):
    '''
    Rename the read IDs in a fastq file with the corresponding sample name.

    Parameters
    fastq : filename or filehandle
        input
    sample_map : dictionary
        entries are {filename: sample}

    yields : SeqRecord
        fastq records
    '''
    try:
        filename_forward, filename_reverse = key
    except ValueError:
        filename_forward = key

    if not reverse:
        direction_flag = 1
        file_to_parse = filename_forward
    else:
        direction_flag = 2
        file_to_parse = filename_reverse

    count = 0
    for record in SeqIO.parse(file_to_parse, 'fastq'):
            desc = '%d#%s/%d' % (count, sample_map[key], direction_flag)
            record.id = desc
            record.description = ''
            count += 1
            yield record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Reformat and concatenate pre-demultiplexed files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('sample', help='sample mapping file')
    parser.add_argument('--output', '-o', default=None, help='name of output file. If there are forward and reverse files, they will be output with a f. and r. prefix respectively. defaults to cat.{sample map file name}.fastq')
    args = parser.parse_args()

    if args.output is None:
        args.output = 'cat.%s%s' % (args.sample, '.fastq')

    out_name = args.output
    forward_out_name = 'f.' + args.output
    reverse_out_name = 'r.' + args.output
    util.check_for_collisions((out_name, forward_out_name, reverse_out_name))

    # parse the sample mapping file
    with open(args.sample, 'r') as f:
        sample_map = sample_file_to_dictionary(f)

    for key in sample_map.keys():
        if isinstance(key, tuple):  # keys of sample_map are tuples when forward and reverse given
            with open(forward_out_name, 'a') as forward_out, open(reverse_out_name, 'a') as reverse_out:
                try:
                    for pair in izip_longest(renamed_fastq_records(key, sample_map), renamed_fastq_records(key, sample_map, reverse=True)):
                        if None in pair:
                            raise ValueError('Files corresponding to the same sample must have the same number of records', key, sample_map[key])
                        SeqIO.write(pair[0], forward_out, 'fastq')
                        SeqIO.write(pair[1], reverse_out, 'fastq')
                except ValueError as e:
                    print '%s. Found unequal number for files %s corresponding to sample %s. Please correct input, delete any output of this script, and try again.' % e.args

        else:  # keys of sample_map strings, only one direction given
            with open(out_name, 'a') as out:
                for record in renamed_fastq_records(key, sample_map):
                    SeqIO.write(record, out, 'fastq')
