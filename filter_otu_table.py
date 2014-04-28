#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from util import *

def parse_args():
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', default = '', help = 'input file')
    parser.add_argument('-o', default = '', help = 'output file')
    parser.add_argument('--row_regex', default = '', help = 'keep rows that match regex')
    parser.add_argument('--col_regex', default = '', help = 'keep cols that match regex')
    parser.add_argument('--transpose', default = False, action = 'store_true', help = 'transpose otu table')
    parser.add_argument('--pseudocount', default = np.nan, type = float, help = 'add pseudocount')
    parser.add_argument('--norm', default = False, action = 'store_true', help = 'normalize')
    parser.add_argument('--log', default = False, action = 'store_true', help = 'log transform')
    parser.add_argument('--locut', default = np.nan, type = float, help = 'lo abundance cutoff (use with max_locut)')
    parser.add_argument('--hicut', default = np.nan, type = float, help = 'hi abundance cutoff (use with max_hicut)')
    parser.add_argument('--max_below_locut', default = np.nan, type = float, help = 'remove otu if (fraction below locut) > max_below_locut')
    parser.add_argument('--min_above_locut', default = np.nan, type = float, help = 'remove otu if (fraction above locut) < min_above_locut')
    parser.add_argument('--max_above_hicut', default = np.nan, type = float, help = 'remove otu if (fraction above hicut) > max_above_hicut')
    parser.add_argument('--min_med', default = np.nan, type = float, help = 'min_med < median < max_med')
    parser.add_argument('--max_med', default = np.nan, type = float, help = 'min_med < median < max_med')
    parser.add_argument('--min_total', default = np.nan, type = float, help = 'remove otu if (total) < min_total')
    parser.add_argument('--top', default = np.nan, type = float, help = 'select most abundant otus (fraction or int)')
    parser.add_argument('--sort', default = False, action = 'store_true', help = 'numeric sort by first column')
    args = parser.parse_args()
    return args


def fmessage(data, text):
    message(text + ', shape = (%d, %d)' %(len(data.index), len(data.columns)))


def filter_otu_table(args, data):
    # filter otu table

    # filter by regex
    if args.row_regex:
        data = data.ix[[bool(re.search(args.row_regex, ri)) for ri in data.index], :]
    if args.col_regex:
        data = data.ix[:, [bool(re.search(args.col_regex, ci)) for ci in data.columns]]
    # transpose
    if args.transpose:
        data = data.transpose()
        fmessage(data, '--transpose: transposing otu table')
    # add pseudocount
    if pd.notnull(args.pseudocount):
        data = data + args.pseudocount
        fmessage(data, '--pseudocount: adding %f to otu table' %(args.pseudocount))
    # normalize
    if args.norm:
        data = data.div(data.sum(axis=1), axis=0)
        fmessage(data, '--norm: normalizing rows of otu table')
    # log transform
    if args.log:
        if pd.isnull(args.pseudocount):
            data = data + 1e-10
        data = np.log(data)
        fmessage(data, '--log: applying log transform')
    # filter by f <= locut
    if pd.notnull(args.locut) and pd.notnull(args.max_below_locut):
        if args.max_below_locut > 1:
            args.max_below_locut = 1.*args.max_below_locut/len(data.index)
        data = data.ix[:, (1.*(data <= args.locut).sum(axis=0) / len(data.index)) < args.max_below_locut]
        fmessage(data, '--locut %f --max_below_locut %f: filtering by minimum abundance' %(args.locut, args.max_below_locut))
    # filter by f > locut
    if pd.notnull(args.locut) and pd.notnull(args.min_above_locut):
        if args.min_above_locut > 1:
            args.min_above_locut = 1.*args.min_above_locut/len(data.index)
        data = data.ix[:, (1.*(data > args.locut).sum(axis=0) / len(data.index)) > args.min_above_locut]
        fmessage(data, '--locut %f --min_above_locut %f: filtering by minimum abundance' %(args.locut, args.min_above_locut))
    # filter by f >= hicut
    if pd.notnull(args.hicut) and pd.notnull(args.max_above_hicut):
        if args.max_above_hicut > 1:
            args.max_above_hicut = 1.*args.max_above_hicut/len(data.index)
        data = data.ix[:, (1.*(data >= args.hicut).sum(axis=0) / len(data.index)) < args.max_above_hicut]
        fmessage(data, '--hicut %f --max_above_hicut %f: filtering by maximum abundance' %(args.hicut, args.max_above_hicut))
    # filter by median
    if pd.notnull(args.min_med):
        data = data.ix[:, data.median(axis=0) >= args.min_med]
        fmessage(data, '--min_med %f: filtering by median abundance' %(args.min_med))
    if pd.notnull(args.max_med):
        data = data.ix[:, data.median(axis=0) <= args.max_med]
        fmessage(data, '--max_med %f: filtering by maximum abundance' %(args.max_med))
    # filter by total
    if pd.notnull(args.min_total):
        data = data.ix[:, data.sum(axis=0) >= args.min_total]
        fmessage(data, '--min_total %f: filtering by total abundance' %(args.min_total))
    # select most abundant otus
    if pd.notnull(args.top):
        if args.top < 1:
            data = data.ix[:, data.median(axis=0).order(ascending=False)[:int(args.top*len(data.index))].index]
            fmessage(data, '--args.top %f: selecting top %f otus' %(args.top))
        elif args.top > 1:
            data = data.ix[:, data.median(axis=0).order(ascending=False)[:int(args.top)].index]
            fmessage(data, '--args.top %d: selecting top %d otus' %(args.top))
    if args.sort == True:
        data.index = data.index.astype(int)
        data = data.sort_index()
    return data


def write_output(args, data):
    # write table as tab-delimited file
    data.to_csv(args.o, sep='\t')


args = parse_args()

if __name__ == '__main__':
    # load input as pandas dataframe
    data = 1.*read_dataframe(args.i)
    fmessage(data, 'loading %s as dataframe' %(args.i))
    data = filter_otu_table(args, data)
    write_output(args, data)
