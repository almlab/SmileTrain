#!/usr/bin/env python

'''
Split an input fastq file into multiple fastq files by
    * opening a cycling set of output file handles
    * iterating over each entry in the input
    * write the entry an output filehandle
    * cycle the output filehandle.

The plus line is stripped of content.

Used for early steps in the pipeline that are embarrassingly parallel.  
'''

import itertools, os, os.path, sys, argparse, itertools, shutil
from Bio import SeqIO
import util

def output_filenames(input_filename, k):
    '''destination filenames foo.fastq.0, etc.'''
    return ['%s.%d' % (input_filename, i) for i in range(k)]

def split_fastq_entries(fastq, filenames):
    '''
    Send entries in the input to filenames, cycling over each filename.
    
    fastq : filename or filehandle
        input
    filenames : list or iterator of strings
        output filenames
        
    returns : nothing
    '''
    
    # open all the filehandles
    fhs = [open(filename, 'w') for filename in filenames]
    fh_cycler = util.cycle(fhs)
    
    # prepare an iterator over the fastq entries
    for record, fh in itertools.izip(SeqIO.parse(fastq, 'fastq'), fh_cycler):
        SeqIO.write(record, fh, 'fastq')
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split a fastq file foo.fastq into multiple fastq files foo.fastq.0, foo.fastq.1, etc.')
    parser.add_argument('fastq', help='input fastq')
    parser.add_argument('n_files', type=int, help='number of split files to output')
    args = parser.parse_args()
    
    filenames = output_filenames(args.fastq, args.n_files)
    util.check_for_collisions(filenames)
    
    if len(filenames) == 1:
        # just copy the file
        shutil.copy(args.fastq, filenames[0])
    else:
        # split the file entry by entry
        split_fastq_entries(args.fastq, filenames)
