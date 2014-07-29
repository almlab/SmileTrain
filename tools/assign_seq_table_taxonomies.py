#!/usr/bin/env python

'''
Get the taxonomies for a seq table using greengenes (or whatever).
'''

import sys, argparse, tempfile, cPickle as pickle
from SmileTrain import ssub

def seq_table_to_fasta(table_fh, fasta_fh):
    '''convert sequence table entries to a fasta'''
    # get the sequences (throwing away the first line "SEQUENCE")
    seqs = [line.split()[0] for line in table_fh]
    seq_header = seqs.pop(0)
    assert(seq_header == "SEQUENCE")
    
    # write the fasta
    for i, seq in enumerate(seqs):
        fasta_fh.write(">seq%s\n%s\n" %(i, seq))
        
def write_tmp_fasta(table_fh, tmp_dir):
    '''make a temporary fasta file'''
    fasta_fh = tempfile.NamedTemporaryFile(suffix='.fst', delete=False, dir=tmp_dir)
    seq_table_to_fasta(table_fh, fasta_fh)
    
    fasta_fn = fasta_fh.name
    fasta_fh.close()
    
    return fasta_fn

def usearch_against_database_cmd(usearch, table_fh, tmp_dir, db, strand='both', sid='0.995'):
    fasta_fn = write_tmp_fasta(table_fh, tmp_dir)
    uc_fh = tempfile.NamedTemporaryFile(suffix='.uc', delete=True, dir=tmp_dir)
    uc_fn = uc_fh.name
    uc_fh.close()
    
    cmd = "%s -usearch_global %s -uc %s -strand %s -id %s -db %s" %(usearch, fasta_fn, uc_fn, strand, sid, db)
    print "  temporary fasta: %s" % fasta_fn
    print "  temporary uc: %s" % uc_fn
    return cmd, uc_fn

def uc_to_ids(uc_fh):
    ids = []
    for i, line in enumerate(uc_fh):
        fields = line.split()
        query = fields[8]
        target = fields[9]
        assert(query == "seq%s" %(i))
        ids.append(target)
    
    return ids

def lookup_taxonomies(ids, tax_pkl_fh):
    d = pickle.load(tax_pkl_fh)
    taxs = [d[qid] for qid in ids]
    return taxs


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('table', help='input sequence table')
    parser.add_argument('db', help='input taxonomy database')
    parser.add_argument('usearch')
    parser.add_argument('tmp_dir')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    
    args = parser.parse_args()
    
    submitter = ssub.Ssub()
    
    cmd, uc_fn = usearch_against_database_cmd(args.usearch, args.table, args.tmp_dir, args.db)
    submitter.submit_and_wait([cmd])
    
    with open(uc_fn) as f:
        ids = uc_to_seq_map(f)
    
    with open(args.db) as f:
        taxs = lookup_taxonomies(ids, f)
        
    args.output.write("\n".join(taxs))