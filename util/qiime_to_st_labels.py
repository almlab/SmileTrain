#!/usr/bin/env python

'''
Convert qiime labels to smile train labels
'''

import argparse, re, sys
from SmileTrain import util

def convert_qiime_label(fasta_entries):
    ''' ">samplename_123 lots of garbage stuff" -> ">samplename;123" '''
    for sid, sequence in fasta_entries:
        # split into "samplename" and "123"
        m = re.match("(.+)_(\d+)", sid)
        if m is None: raise RuntimeError("line %s did not parse" %(sid))
        sample_name, sequence_number = m.groups()
        
        # create a new sid
        new_sid = "sample=%s;%s" %(sample_name, sequence_number)
        
        yield new_sid, sequence
    

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input fasta')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output fasta (default stdout)')
    args = parser.parse_args()
    
    with open(args.input) as f:
        for sid, seq in convert_qiime_label(util.fasta_entries(f)):
            args.output.write(">%s\n%s\n" %(sid, seq))
        