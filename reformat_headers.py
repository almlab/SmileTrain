#!/usr/bin/env python

'''
Reformat headers when user started with already demultiplexed data

    @OURSEQ:lolapalooza1234#MySample/1
    AACCGGTT
    +
    abcdefgh

becomes output like
    @sample=MySample;1
    AACCGGTT
    +whatever
    abcdefgh
    
where the ;1 means it's the first read that mapped to donor1_day5.
'''

import sys, argparse, string, itertools, re
from Bio import SeqIO

def renamed_fastq_records(fastq):
    '''
    Rename the read IDs in a fastq file with the corresponding sample name. Get the sample name
    from the header itself

    Parameters
    fastq : filename or filehandle
        input

    yields : SeqRecord
        fastq records
    '''

    sample_counts = {}

    for record in SeqIO.parse(fastq, 'fastq'):
        # look for the barcode from the read ID line
        m = re.match('.*#(.+)\/(\d)$', record.id)
        sample = m.group(1)
        read_direction = m.group(2)
        
        if sample in sample_counts:
            sample_counts[sample] += 1
        else:
            sample_counts[sample] = 1

        record.id = "sample=%s;%d/%s" %(sample, sample_counts[sample], read_direction)
        
        # expunge other parts of title
        record.description = ''
        yield record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Demultiplex fastq entries by barcode', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fastq', help='input fastq file')  
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output fastq')
    args = parser.parse_args()

    # get a set of reads
    for record in renamed_fastq_records(args.fastq):
        SeqIO.write(record, args.output, 'fastq')