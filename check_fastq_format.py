#!/usr/bin/env python

'''
Before running any steps in the pipeline, check your data to see if the fastq input is in
Illumina 1.3-1.7 format, specificially that
    * the @ line follows the format @stuff:with:colons#BARCODE/1 (or /2 for reverse reads)
    * the quality line uses Sanger encoding (ascii offset 64, high score h to low score B)
'''

import itertools, os, os.path, sys, argparse, shutil, re
from Bio import SeqIO
import util

def parse_fastq_record_id(record):
    '''BioPython fastq record "@lol/1" -> "lol"'''
    m = re.match('^(.+)/[12]', record.id)
    if m is None:
        raise RuntimeError("fastq record line did not parse: %s" % record.id)
    else:
        rid = m.group(1)
        
    return rid

def check_illumina_format(inputs, targets):
    '''
    Raise error if any of the input filenames are not in the desired format.
    
    inputs : filenames or filehandle, or iterable of strings/filehandles
        files to be checked
    target : string, or iterable of strings
        'illumina13' or 'illumina18' or 'ambiguous', or some combination
    '''
    
    # convert to lists if required
    inputs = util.listify(inputs)
    targets = util.listify(targets)
    
    # check that input targets are in acceptable list
    acceptable_targets = set(['illumina13', 'illumina18', 'ambiguous'])
    if not set(targets).issubset(acceptable_targets):
        bad_targets = targets - acceptable_targets
        raise ArgumentError("unrecognized format type(s): %s" % bad_targets)
    
    # check all the formats
    formats = [file_format(i) for i in inputs]
    tests = [form in targets for form in formats]
    bad_files = [i for i, test in zip(i, tests) if test == False]
    bad_forms = [form for form, test in zip(formats, tests) if test == False]
    
    # complain if something went wrong
    if False in tests:
        bad_info = "\n".join([" ".join([i, form]) for i, form in zip(bad_files, bad_forms)])
        raise RuntimeError("files do not appear to be in %s format: \n%s" % (targets, bad_info))
    else:
        return True

def file_format(fastq, max_entries=10):
    '''
    what format is the file?
    
    returns : string
        'illumina13', 'illumin18', or 'ambiguous'
    '''
    
    for i, record in enumerate(SeqIO.parse(fastq, 'fastq')):
        if i > max_entries:
            raise RuntimeError("could not verify format after %d entries" % max_entries)
        
        # make sure we can parse the at line
        rid = parse_fastq_record_id(record)
        
        # check the quality line's character content
        return fastq_record_format(record)
    
    raise RuntimeError("fell off end")
        
def fastq_record_format(record):
    '''
    Guess the fastq format using the quality codes.
    
    record : BioPython Seq object
        to be analyzed
        
    returns : 'illumina13', 'illumina18', 'ambiguous'
        either Illumina 1.3-1.7 (B--h), Illumina 1.8 ("--J), or ambiguous
    '''
    
    scores = record.letter_annotations['phred_quality']
    min_score = min(scores)
    max_score = max(scores)
    
    if 1 <= min_score and min_score < 32 and max_score <= 41:
        return 'illumina18'
    elif 32 <= min_score and 41 < max_score and max_score <= 71:
        return 'illumina13'
    elif min_score < 1 or max_score > 71:
        raise RuntimeError("quality scores don't match known encoding: min=%s max=%s" %(min_score, max_score))
    else:
        return 'ambiguous'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verify that fastq is in Illumina 1.3-1.7 format')
    parser.add_argument('fastq', help='input fastq')
    args = parser.parse_args()
    
    with open(args.fastq) as f:
        format_guess = file_format(f)
        if format_guess == 'illumina13':
            print "Looks like Illumina 1.3-1.7 format. Proceed with the pipeline!"
        elif format_guess == 'illumina18':
            print "Looks like Illumina 1.8 format. You may need to convert. Beware..."
        elif format_guess == 'ambiguous':
            print "Could be either 1.3-1.7 or 1.8 format. Ambiguous. Proceed with caution."
        else:
            raise RuntimeError