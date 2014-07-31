#!/usr/bin/env python

'''
Get greengenes taxonomies. Given an otu table with otu ids in the first column, search through
the greengenes taxonomy list. Output the taxonomies in order.

If the input database is a pickle, just load that dictionary.
'''

import sys, argparse, re, cPickle as pickle

def table_ids(fn):
    '''get the otu ids from the otu table with filename fn'''
    with open(fn) as f:
        ids = [line.split()[0] for line in f]

    # remove the first item, which is "OTU_ID"
    ids.pop(0)

    return ids

def uc_ids(fn):
    with open(fn) as f:
        ids = [line.split()[-1] for line in f]
        
    return ids

def list_ids(fn):
    with open(fn) as f:
        ids = [line.strip() for line in f]
    
    return ids

def taxa_dictionary(fn, ids):
    '''get the second field in lines whose first fields match ids'''
    # populate a hash otu_id => taxonomy
    d = {}
    with open(fn) as f:
        for line in f:
            fields = line.split()
            otu = fields[0]
            tax = " ".join(fields[1:])
            if otu in ids:
                d[otu] = tax

    return d


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('db', help='input database')
    input_type = parser.add_mutually_exclusive_group(required=True)
    input_type.add_argument('-t', '--table', default=None, help='input otu table')
    input_type.add_argument('-u', '--uc', default=None, help='input uc file')
    input_type.add_argument('-l', '--list', default=None, help='input plain text list')
    
    parser.add_argument('-i', '--no_match_id', default=None, help='OTU ID for no match (default "no_match" for table/list; "*" for uc)')
    parser.add_argument('-x', '--no_match_tax', default='k__; p__; c__; o__; f__; g__; s__', help='taxonomy for unmatched OTU ID (default is QIIME taxonomy format)')
    args = parser.parse_args()
    
    # depending on the input type, adjust the no match id and parsing function
    if args.table is not None:
        ids = table_ids(args.table)
        if args.no_match_id is None: args.no_match_id = 'no_match'
    elif args.uc is not None:
        ids = uc_ids(args.uc)
        if args.no_match_id is None: args.no_match_id = '*'
    elif args.list is not None:
        ids = list_ids(args.list)
        if args.no_match_id is None: args.no_match_id = 'no_match'

    # check if the database file ends in .pkl or .pickle
    # if it is, used a pickled dictionary
    # otherwise, just search line by line
    if re.search("\.(pkl|pickle)$", args.db):
        with open(args.db, 'rb') as f:
            d = pickle.load(f)
    else:   
        d = taxa_dictionary(args.db, ids)
        
    d[args.no_match_id] = args.no_match_tax
        
    print "\n".join([d[i] for i in ids])