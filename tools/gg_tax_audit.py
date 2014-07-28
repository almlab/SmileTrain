#!/usr/bin/env python

'''
A pipeline for examining the taxonomies of some particular OTUs from greengenes
'''

import argparse, os, ConfigParser
from SmileTrain import ssub, util


def matching_labels(uc, target):
    '''get the seq ids that match some greengenes id'''
    labels = []
    for line in uc:
        fields = line.split()
        label = fields[8]
        this_target = fields[9]
        
        if this_target == target:
            labels.append(label)
            
    return labels

def matching_entries(fasta, labels):
    '''get the sequences that have some label'''
    entries = [entry for entry in util.fasta_entries(fasta) if entry[0] in labels]
    return entries
    

if __name__ == '__main__':
    # open the config file sister to this script
    config = ConfigParser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), os.pardir, 'user.cfg'))
    
    parser = argparse.ArgumentParser()
    parser.add_argument('gg_id', help='reference id')
    parser.add_argument('uc', help='input uc file')
    parser.add_argument('gg_rep_set', help='input gg rep set fasta')
    parser.add_argument('--fasta', default='q.derep.fst', help='dereplicated fasta (default: q.derep.fst)')
    parser.add_argument('top_match', help='top match fasta')
    parser.add_argument('matches', help='output matched seqs fasta')
    parser.add_argument('-m', '--mismatches', default=None, help='output mismatch file (default: none)')
    parser.add_argument('-u', '--out_uc', default='out.uc', help='output uc (default: out.uc)')
    parser.add_argument('--strand', default='both', help='strand (default: both)')
    parser.add_argument('--sid', default='0.995', help='sequence identity (default: 0.995)')
    parser.add_argument('--tax', default='tax.txt', help='taxonomy output (default: tax.txt)')
    parser.add_argument('--dry', action='store_true', help='don\'t submit jobs, just show them?')
    
    args = parser.parse_args()
    
    # get some values from the config file
    usearch = config.get('User', 'usearch')
    gg_tax = config.get('Data', 'greengenes_taxonomy')
    library = config.get('Scripts', 'library')
    
    # get the matching seq labels from the uc file
    util.message('getting matching seq labels from uc...')
    with open(args.uc) as f:
        labels = matching_labels(f, args.gg_id)
    
    # get the sequences matching these labels
    util.message('getting matching sequences from fasta...')
    with open(args.fasta) as f:
        entries = matching_entries(f, labels)
        
    # write out the matches
    util.message('writing matches...')
    with open(args.matches, 'w') as f:
        f.write(util.fasta_entries_to_string(entries))
        
    with open(args.top_match, 'w') as f:
        f.write(util.fasta_entry_list_to_string(entries[0]) + "\n")
    
    # start an internal ssub object
    submitter = ssub.Ssub()
    
    # count the mismatches if desired
    if args.mismatches is not None:
        util.message('counting mismatches...')
        cmd = 'python %s/tools/count_mismatches.py %s > %s' %(library, args.matches, args.mismatches)
        submitter.submit([cmd])
    
    cmd1 = '%s -search_global %s -db %s -uc %s -uc_allhits -maxaccepts 0 -maxrejects 0 -strand %s -id %s' %(usearch, args.top_match, args.gg_rep_set, args.out_uc, args.strand, args.sid)
    cmd2 = '%s/tools/get_taxonomies.py %s --uc %s > %s' %(library, gg_tax, args.out_uc, args.tax)
    if args.dry:
        print "run these commands on a node:"
        print "%s; %s" %(cmd1, cmd2)
    else:
        util.message('searching for sequences...')
        submitter.submit_and_wait([cmd1])

        util.message('matching taxonomies...')
        submitter.submit_and_wait([cmd])