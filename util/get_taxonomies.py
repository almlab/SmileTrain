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

def matching_fields(fn, ids):
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

    # get the taxonomies in the order that the ids were given
    taxa = [d[i] for i in ids]
    return taxa

def pickled_taxa(fn, ids):
    '''get matching taxa from a pickled tax dictionary'''
    with open(fn, 'rb') as f:
        d = pickle.load(f)
    
    taxa = [d[i] for i in ids]
    return taxa


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input otu table')
    parser.add_argument('db', help='input database')
    args = parser.parse_args()

    ids = table_ids(args.input)

    # check if the database file ends in .pkl or .pickle
    # if it is, used a pickled dictionary
    # otherwise, just search line by line
    if re.search("\.(pkl|pickle)$", args.db):
        taxa = pickled_taxa(args.db, ids)
    else:   
        taxa = matching_fields(args.db, ids)
        
    print "\n".join(taxa)
