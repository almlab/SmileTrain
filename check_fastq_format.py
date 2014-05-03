#!/usr/bin/env python

'''
Before running any steps in the pipeline, check your data to see if the fastq input is in
Illumina 1.3-1.7 format, specificially that
    * the @ line follows the format @stuff:with:colons#BARCODE/1 (or /2 for reverse reads)
    * the quality line uses Sanger encoding (ascii offset 64, high score h to low score B)
'''

import itertools, os, os.path, sys, argparse, shutil
import util

def check_illumina_format(filenames, target):
    '''
    Raise error if any of the input filenames are not in the desired format.
    
    filenames : iterable of strings
        files to be checked
    target : string
        'illumina13' or 'illumina18'
    '''
    
    if target not in ['illumina13', 'illumina18']:
        raise ArgumentError("unrecognized format type: %s" % target)
    
    format_f = lambda fn: file_format(open(fn))
    
    formats = [file_format(open(fn)) for fn in filenames]
    tests = [form == target for form in formats]
    bad_files = [fn for fn, test in zip(filenames, tests) if test == False]
    bad_forms = set([form for form, test in zip(formats, tests) if test == False])
    
    if False in tests:
        bad_info = "\n".join([" ".join([fn, form]) for fn, form in zip(filenames, bad_forms)])
        raise RuntimeError("files do not appear to be in %s format: \n%s" % (target, bad_info))

def file_format(fastq, max_entries=10):
    '''
    what format is the file?
    
    returns : string
        'illumina13', 'illumin18', or 'ambiguous'
    '''
    
    for i, (at_line, seq_line, qua_line) in enumerate(util.fastq_iterator(fastq)):
        if i > max_entries:
            raise RuntimeError("could not verify format after %d entries" % max_entries)
        
        # make sure we can parse the at line
        rid = util.fastq_at_line_to_id(at_line)
        
        # check the quality line's character content
        return quality_line_format(qua_line)
        
def quality_line_format(line):
    '''
    Guess the fastq format using the quality codes.
    
    line : string
        quality line from a fastq entry
        
    returns : 'illumina13', 'illumina18', 'ambiguous'
        either Illumina 1.3-1.7 (B--h), Illumina 1.8 ("--J), or ambiguous
    '''
    
    matches_illumina13 = all([c in util.illumina13_codes for c in list(line)])
    matches_illumina18 = all([c in util.illumina18_codes for c in list(line)])
    
    if matches_illumina13 and matches_illumina18:
        return 'ambiguous'
    elif matches_illumina13 and not matches_illumina18:
        return 'illumina13'
    elif matches_illumina18 and not matches_illumina13:
        return 'illumina18'
    elif not matches_illumina18 and not matches_illumina13:
        raise RuntimeError("quality line doesn't match known encoding: %s" % line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verify that fastq is in Illumina 1.3-1.7 format')
    parser.add_argument('fastq', help='input fastq')
    args = parser.parse_args()
    
    with open(args.fastq) as f:
        if is_illumina13_format(f):
            print "Looks like Illumina 1.3-1.7 format. Proceed with the pipeline!"
        else:
            print "Doesn't look like Illumina 1.3-1.7 format. Investigate before proceeding."