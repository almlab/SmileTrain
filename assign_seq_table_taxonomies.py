#!/usr/bin/env python

'''
Get the taxonomies for a seq table using greengenes (or whatever).
'''

import sys, argparse, tempfile, cPickle as pickle, ConfigParser, os
from Bio import Seq, SeqIO, SeqRecord
from SmileTrain import ssub

def seq_table_to_fasta(table_fh, fasta_fh):
    '''convert sequence table entries to a fasta'''
    # get the sequences (throwing away the first line "SEQUENCE")
    seqs = [line.split()[0] for line in table_fh]
    seq_header = seqs.pop(0)
    assert(seq_header.lower() == "sequence")
    
    # write the fasta
    for i, seq in enumerate(seqs):
        record = SeqRecord.SeqRecord(Seq.Seq(seq), id="seq%s" %(i))
        SeqIO.write(record, fasta_fh, 'fasta')
        
def write_tmp_fasta(table_fh, tmp_dir):
    '''make a temporary fasta file'''
    fasta_fh = tempfile.NamedTemporaryFile(suffix='.fst', delete=False, dir=tmp_dir)
    seq_table_to_fasta(table_fh, fasta_fh)
    
    fasta_fn = fasta_fh.name
    fasta_fh.close()
    
    return fasta_fn

def usearch_against_database_cmd(usearch, table_fh, tmp_dir, db, fid, no_hit=None, strand='both'):
    fasta_fn = write_tmp_fasta(table_fh, tmp_dir)
    uc_fh = tempfile.NamedTemporaryFile(suffix='.uc', delete=True, dir=tmp_dir)
    uc_fn = uc_fh.name
    uc_fh.close()
    
    cmd = "%s -usearch_global %s -uc %s -strand %s -id %s -db %s" %(usearch, fasta_fn, uc_fn, strand, fid, db)
    
    if no_hit is not None:
        cmd += " -notmatched %s" %(no_hit)
    
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

def lookup_taxonomies(ids, tax_pkl_fh, no_match=None, no_match_label=None):
    d = pickle.load(tax_pkl_fh)
    
    if no_match is not None:
        if no_match_label is None:
            raise RuntimeError("no match symbol specified (%s) but not label!" %(no_match))
        
        taxs[no_match] = no_match_label
    
    taxs = [d[qid] for qid in ids]
    return taxs


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('table', help='input sequence table')
    parser.add_argument('-s', '--sid', default='99', help='greengenes repset identity (default: 99)')
    parser.add_argument('-i', '--fid', default='0.995', help='fractional identity to make hit (default: 0.995)')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output file (default stdout)')
    parser.add_argument('--no_hit', default=None, help='get unmatched sequences as separate fasta?')
    parser.add_argument('--no_match', '-n', default='*', help='no match indicator in uc (default: *)')
    parser.add_argument('--no_match_label', '-l', default='k__; p__; c__; o__; f__; g__; s__', help='taxonomy for no match (default: Qiime empty)')
    
    args = parser.parse_args()
    
    # grab the configuration file
    config = ConfigParser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'user.cfg'))
    
    usearch = config.get('User', 'usearch')
    gg_dir = config.get('Data', 'greengenes')
    gg_tax = config.get('Data', 'greengenes_taxonomy')
    library = config.get('Scripts', 'library')
    tmp_dir = config.get('User', 'tmp_directory')
    
    submitter = ssub.Ssub()
    
    gg_fasta = os.path.join(gg_dir, "%s_otus.fasta" %(args.sid))
    with open(args.table) as f:
        cmd, uc_fn = usearch_against_database_cmd(usearch, f, tmp_dir, gg_fasta, args.fid, args.no_hit)
    submitter.submit_and_wait([cmd])
    
    with open(uc_fn) as f:
        ids = uc_to_ids(f)
    
    with open(gg_tax) as f:
        taxs = lookup_taxonomies(ids, f, args.no_match, args.no_match_label)
        
    args.output.write("\n".join(taxs) + "\n")