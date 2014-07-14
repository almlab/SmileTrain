#!/usr/bin/env python

'''
Pickle greengenes taxonomies. This makes taxonomy lookups fast & easy.
'''

import argparse, cPickle as pickle

def taxonomy_file_to_dict(fn):
    '''tab-separated taxonomy table to dictionary'''
    d = {}
    with open(fn) as f:
        for line in f: 
            otu, tax = line.strip().split("\t")
            d[otu] = tax
    
    return d


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input taxonomy table')
    parser.add_argument('-o', '--output', default='taxonomy.pkl', help='output pickle file (default: taxonomy.pkl')
    args = parser.parse_args()
    
    tax_dict = taxonomy_file_to_dict(args.input)
    
    # write the dictionary as a binary pickle
    with open(args.output, 'wb') as f:
        pickle.dump(tax_dict, f)