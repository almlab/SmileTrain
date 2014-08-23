#!/usr/bin/env python

'''
Main script in the pipeline. Produces lists of commands and submits them using ssub.
Options allow the user to run individual parts of the pipeline or the entire thing.
'''

import argparse, os, ConfigParser
import ssub
import util, check_fastq_format
from util import *

# open the config file sister to this script
config = ConfigParser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'user.cfg'))

def parse_args():
    '''arguments to be parsed and passed to the OTU_Caller object'''

    # create argument parser
    parser = argparse.ArgumentParser(description=util.banner, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # add groups
    group1 = parser.add_argument_group('Pipeline')
    group2 = parser.add_argument_group('Input files')
    group3 = parser.add_argument_group('Convert format')
    group4 = parser.add_argument_group('Remove primers')
    group5 = parser.add_argument_group('Merge reads')
    group6 = parser.add_argument_group('Demultiplex')
    group7 = parser.add_argument_group('Quality filtering')
    group8 = parser.add_argument_group('Dereplicate')
    group9 = parser.add_argument_group('Chimeras')
    group10 = parser.add_argument_group('Indexing')
    group11 = parser.add_argument_group('Clustering')
    group12 = parser.add_argument_group('Options')
    
    # add arguments
    group1.add_argument('--all', action='store_true', help='Run primer, merge, demultiplex, filter, derep, index, ref_gg, and otu?')
    group1.add_argument('--check', action='store_true', help='Check input file format?')
    group1.add_argument('--split', action='store_true', help='Split the fastq files?')
    group1.add_argument('--convert', action='store_true', help='Convert fastq format?')
    group1.add_argument('--primers', action='store_true', help='Remove primers?')
    group1.add_argument('--merge', action='store_true', help='Merge forward and reverse reads?')
    group1.add_argument('--demultiplex', default = False, action = 'store_true', help = 'Demultiplex?')
    group1.add_argument('--qfilter', default = False, action = 'store_true', help = 'Quality filter?')
    group1.add_argument('--chimeras', default = False, action = 'store_true', help = 'Chimera slay?')
    group1.add_argument('--dereplicate', action='store_true', help='Dereplicate?')
    group1.add_argument('--index', action='store_true', help='Make index file?')
    group1.add_argument('--denovo', default = False, action = 'store_true', help = 'Denovo clustering (UPARSE)?')
    group1.add_argument('--ref_gg', default = False, action = 'store_true', help = 'Reference mapping (Greengenes)?')
    group1.add_argument('--open_ref_gg', action='store_true', help='Reference map (Greengenes) and then denovo cluster?')
    group1.add_argument('--seq_table', action='store_true', help='Make a sequence table?')
    group1.add_argument('--seq_tax', action='store_true', help='Get taxonomies for sequence table?')
    group1.add_argument('--otu_table', action='store_true', help='Make OTU table?')
    group2.add_argument('--forward', '-f', help='Input fastq (forward)')
    group2.add_argument('--reverse', '-r', help='Input fastq (reverse)')
    group2.add_argument('-p', help='Primer sequence (forward)')
    group2.add_argument('-q', help='Primer sequence (reverse)')
    group2.add_argument('--barcodes', '-b', default=None, help='Barcodes list')
    group4.add_argument('--p_mismatch', default=1, type=int, help='Number of mismatches allowed in primers')
    group6.add_argument('--b_mismatch', default=1, type=int, help='Number of mismatches allowed in barcodes')
    group6.add_argument('--merged', action='store_true', help='Files were merged in a previous step?')
    group7.add_argument('--truncqual', default = 2, type = int, help = '')
    group7.add_argument('--maxee', default = 2., type = float, help = 'Maximum expected error (UPARSE)')
    group9.add_argument('--gold_db', default=config.get('Data', 'gold'), help='Gold 16S database')
    group11.add_argument('--sids', default='91,94,97,99', help='Sequence identities for clustering')
    group12.add_argument('--n_cpus', '-n', default = 1, type = int, help='Number of CPUs')
    group12.add_argument('--dry_run', '-z', action='store_true', help='submit no jobs; suppress file checks; ust print output commands')
    
    # parse arguments
    if __name__ == '__main__':
        args = parser.parse_args()
    else:
        args = parser.parse_args('')
    
    # process arguments
    if args.all == True:
        args.check = args.split = args.convert = args.primers = args.merge = args.demultiplex = args.qfilter = args.dereplicate = args.index = args.ref_gg = args.otu_table = True
    args.sids = map(int, args.sids.split(','))
    
    # check combinations
    if args.demultiplex and args.barcodes is None:
        raise RuntimeError("--demultiplex selected but no barcode mapping file specified")
    
    if args.check or args.split or args.convert or args.primers or args.merge or args.demultiplex or args.qfilter:
        if args.forward is None and args.reverse is None:
            raise RuntimeError("no fastq files selected")
        
    return args


class OTU_Caller():
    '''
    A namespace and method container for the parsed command line options. Submitted jobs
    refer to other scripts in the library.
    '''

    def __init__(self):
        # initialize variables
        self.usearch = config.get('User', 'usearch')
        self.ggdb = config.get('Data', 'greengenes')
        self.library = config.get('Scripts', 'library')

        # copy command line arguments
        self.__dict__.update(parse_args().__dict__)
        # create filenames
        self.get_filenames()

        # create an internal ssub object. this way we can set ssub's options without
        # relying on the command line.
        self.ssub = ssub.Ssub()
        self.ssub.n_cpus = self.n_cpus
    
    def get_filenames(self):
        '''Generate filenames to use in pipeline'''

        if self.forward:
            f_base = os.path.basename(self.forward)
        if self.reverse:
            r_base = os.path.basename(self.reverse)
            
        self.fi = ['%s.%d' %(self.forward, i) for i in range(self.n_cpus)] # forward reads (split)
        self.ri = ['%s.%d' %(self.reverse, i) for i in range(self.n_cpus)] # reverse reads (split)
        self.mi = ['%s.%d.merge' %(self.forward, i) for i in range(self.n_cpus)] # merged reads (split)
        self.Fi = ['%s.%d.tmp' %(self.forward, i) for i in range(self.n_cpus)] # forward reads (temp)
        self.Ri = ['%s.%d.tmp' %(self.reverse, i) for i in range(self.n_cpus)] # reverse reads (temp)
        self.Mi = ['%s.%d.tmp' %(self.forward, i) for i in range(self.n_cpus)] # merged reads (temp)
        self.ci = ['q.%d.fst' %(i) for i in range(self.n_cpus)] # current reads
        self.Ci = ['q.%d.tmp' %(i) for i in range(self.n_cpus)] # current reads (temp)
        self.oi = ['otus.%d.fst' %(sid) for sid in self.sids] # otu representative sequences
        self.Oi = ['otus.%d.tmp' %(sid) for sid in self.sids] # otu representative sequences (temp)
        self.uc = ['otus.%d.uc' %(sid) for sid in self.sids] # uclust output files
        self.xi = ['otus.%d.counts' %(sid) for sid in self.sids] # otu tables (counts)
        
        self.open_fst = ['q.%d.no_match.fst' %(sid) for sid in self.sids] # otu rep seqs
        self.orig_uc = ['otus.%d.ref_db.uc' %(sid) for sid in self.sids] # uclust output for initial reference-based mapping
        self.open_uc = ['otus.%d.no_match.uc' %(sid) for sid in self.sids] # uclust output files, for open reference clustering
        self.open_otu_tmp = ['otus.%d.no_match.tmp' %(sid) for sid in self.sids] # otu rep seqs (tmp)
        self.open_otu = ['otus.%d.no_match.fst' %(sid) for sid in self.sids] # otu rep seqs
        self.seq_tax_fn = 'seq.tax'
        
        # if reference-based clustering at 97%, make sure reads are 98.5% dissimilar
        self.reference_map_sids = [0.5*(100.0+float(sid)) for sid in self.sids]
        self.reference_map_pcts = [100.0 - sid for sid in self.reference_map_sids]
        
        # Get database for read mapping
        if self.denovo == True:
            self.db = self.oi
        elif self.ref_gg == True:
            self.db = ['%s/%d_otus.fasta' %(self.ggdb, sid) for sid in self.sids]
            
    def check_format(self):
        '''Make sure we have the correct input format'''
        files = []
        if self.forward:
            files.append(self.forward)
        
        if self.reverse:
            files.append(self.reverse)
            
        message('Testing format of %s' %(" ".join(files)))
        
        if self.dry_run:
            cmds = ['python %s/check_fastq_format.py %s' %(self.library, f) for f in files]
            message("\n".join(cmds), indent=0)
        else:
            check_fastq_format.check_illumina_format(files, 'illumina13')
    
    def split_fastq(self):
        '''Split forward and reverse reads (for parallel processing)'''

        # do forward only if there is a forward file; similar for reverse
        do_forward = self.forward
        do_reverse = self.reverse

        # check for inputs and collisions of output
        if do_forward:
            util.check_for_nonempty(self.forward, self.dry_run)
            util.check_for_collisions(['%s.%s' %(self.forward, i) for i in range(self.n_cpus)], self.dry_run)
        if do_reverse:
            util.check_for_nonempty(self.reverse, self.dry_run)
            util.check_for_collisions(['%s.%s' %(self.reverse, i) for i in range(self.n_cpus)], self.dry_run)
        
        # Get list of commands
        cmds = []
        if do_forward:
            cmd = 'python %s/split_fastq.py %s %s' %(self.library, self.forward, self.n_cpus)
            cmds.append(cmd)
        if do_reverse:
            cmd = 'python %s/split_fastq.py %s %s' %(self.library, self.reverse, self.n_cpus)
            cmds.append(cmd)
        
        # submit commands
        self.ssub.submit_and_wait(cmds, out=self.dry_run)

        # validate output
        if do_forward:
            util.check_for_nonempty(self.fi, self.dry_run)
        if do_reverse:
            util.check_for_nonempty(self.ri, self.dry_run)
            
    def convert_format(self):
        '''Convert to compatible fastq format'''
        
        if self.forward:
            util.check_for_nonempty(self.fi, self.dry_run)
            util.check_for_collisions(self.Fi, self.dry_run)
            
        if self.reverse:
            util.check_for_nonempty(self.ri, self.dry_run)
            util.check_for_collisions(self.Fi, self.dry_run)
            
        cmds = []
        for i in range(self.n_cpus):
            if self.forward:
                cmd = 'python %s/convert_fastq.py %s > %s' %(self.library, self.fi[i], self.Fi[i])
                cmds.append(cmd)
            if self.reverse:
                cmd = 'python %s/convert_fastq.py %s > %s' %(self.library, self.ri[i], self.Ri[i])
                cmds.append(cmd)
                
        self.ssub.submit_and_wait(cmds, out=self.dry_run)
        
        # validate output and move files
        if self.forward:
            util.check_for_nonempty(self.Fi, dry_run=self.dry_run)
            self.ssub.move_files(self.Fi, self.fi, out=self.dry_run)
            util.check_for_nonempty(self.fi, dry_run=self.dry_run)

        if self.reverse:
            util.check_for_nonempty(self.Ri, dry_run=self.dry_run)
            self.ssub.move_files(self.Ri, self.ri, out=self.dry_run)
            util.check_for_nonempty(self.ri, dry_run=self.dry_run)
    
    def remove_primers(self):
        '''Remove diversity region + primer and discard reads with > 2 mismatches'''

        # do forward only if there is a forward read file and a forward primer
        # similar for reverse
        do_forward = self.forward and self.p
        do_reverse = self.reverse and self.q

        # check that something is being done
        if not (do_forward or do_reverse):
            raise RuntimeError("remove primers called with bad input: a file or primer is missing")

        # check for inputs and collisions of output
        if do_forward:
            util.check_for_nonempty(self.fi, self.dry_run)
            util.check_for_collisions(self.Fi, self.dry_run)
        if do_reverse:
            util.check_for_nonempty(self.ri, self.dry_run)
            util.check_for_collisions(self.Ri, self.dry_run)
        
        # get list of commands using forward, reverse, or both
        cmds = []
        for i in range(self.n_cpus):
            if do_forward:
                cmd = 'python %s/remove_primers.py %s %s --max_primer_diffs %d > %s' %(self.library, self.fi[i], self.p, self.p_mismatch, self.Fi[i])
                cmds.append(cmd)
            if do_reverse:
                cmd = 'python %s/remove_primers.py %s %s --max_primer_diffs %d > %s' %(self.library, self.ri[i], self.q, self.p_mismatch, self.Ri[i])
                cmds.append(cmd)
        
        # submit commands
        self.ssub.submit_and_wait(cmds, out=self.dry_run)

        # validate output and move files
        if do_forward:
            util.check_for_nonempty(self.Fi, dry_run=self.dry_run)
            self.ssub.move_files(self.Fi, self.fi, out=self.dry_run)
            util.check_for_nonempty(self.fi, dry_run=self.dry_run)

        if do_reverse:
            util.check_for_nonempty(self.Ri, dry_run=self.dry_run)
            self.ssub.move_files(self.Ri, self.ri, out=self.dry_run)
            util.check_for_nonempty(self.ri, dry_run=self.dry_run)
    
    def merge_reads(self):
        '''Merge forward and reverse reads using USEARCH'''

        # check for inputs and collisions
        util.check_for_nonempty(self.fi + self.ri, dry_run=self.dry_run)
        util.check_for_collisions(self.Fi + self.Ri, dry_run=self.dry_run)
        
        # check that usearch is ready to go
        assert(util.is_executable(self.usearch))
        
        # Intersect forward and reverse reads
        cmds = []
        for i in range(self.n_cpus):
            cmd = 'python %s/intersect.py %s %s %s %s' %(self.library, self.fi[i], self.ri[i], self.Fi[i], self.Ri[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, out = self.dry_run)
        
        # make sure there was a nonempty result. copy to new location, and make sure the copy worked
        util.check_for_nonempty(self.Fi + self.Ri, dry_run=self.dry_run)
        self.ssub.move_files(self.Fi + self.Ri, self.fi + self.ri, out = self.dry_run)
        util.check_for_nonempty(self.fi + self.ri, dry_run=self.dry_run)
        
        # Merge reads
        cmds = []
        for i in range(self.n_cpus):
            cmd = '%s -fastq_mergepairs %s -reverse %s -fastq_truncqual %d -fastqout %s' %(self.usearch, self.fi[i], self.ri[i], self.truncqual, self.mi[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, out = self.dry_run)
        
        util.check_for_nonempty(self.mi, dry_run=self.dry_run)
        self.ssub.remove_files(self.fi + self.ri, out = self.dry_run)
    
    def demultiplex_reads(self):
        '''Demultiplex samples using index and barcodes'''
        
        util.check_for_nonempty(self.ci, dry_run=self.dry_run)
        util.check_for_collisions(self.Ci, dry_run=self.dry_run)

        cmds = []
        for i in range(self.n_cpus):
            cmd = 'python %s/map_barcodes.py %s %s --max_barcode_diffs %d > %s' %(self.library, self.ci[i], self.barcodes, self.b_mismatch, self.Ci[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        
        util.check_for_nonempty(self.Ci, dry_run=self.dry_run)
        self.ssub.move_files(self.Ci, self.ci, self.dry_run)
        util.check_for_nonempty(self.ci, dry_run=self.dry_run)
    
    def quality_filter(self):
        '''Quality filter with truncqual and maximum expected error'''
        
        # validate input/output
        util.check_for_nonempty(self.ci, self.dry_run)
        util.check_for_collisions(self.Ci, self.dry_run)
        
        # check that usearch is ready to go
        assert(util.is_executable(self.usearch))
        
        # check that the files are in the right format
        if self.dry_run:
            cmds = ['python %s/check_fastq_format.py %s' %(self.library, f) for f in self.ci]
            message("\n".join(cmds), indent=0)
        else:
            check_fastq_format.check_illumina_format(self.ci, ['illumina18', 'ambiguous'])

        cmds = []
        for i in range(self.n_cpus):
            cmd = '%s -fastq_filter %s -fastq_truncqual %d -fastq_maxee %f -fastaout %s' %(self.usearch, self.ci[i], self.truncqual, self.maxee, self.Ci[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        
        cmd = 'cat %s > q.fst' %(' '.join(self.Ci))
        self.ssub.run_local([cmd], out = self.dry_run)
        util.check_for_nonempty('q.fst', dry_run=self.dry_run)
        
        cmd = 'rm %s' %(' '.join(self.Ci))
        self.ssub.run_local([cmd], out = self.dry_run)
    
    def dereplicate_reads(self):
        '''Concatenate files and dereplicate'''

        cmd = 'python %s/derep_fulllength.py q.fst -o q.derep.fst' %(self.library)
        self.ssub.submit_and_wait([cmd], self.dry_run)
        util.check_for_nonempty('q.derep.fst', dry_run=self.dry_run)
        
    def make_index(self):
        '''Make an index file'''
        
        # verify input & check for collisions with output
        util.check_for_nonempty(['q.fst', 'q.derep.fst'], dry_run=self.dry_run)
        util.check_for_collisions('q.index', self.dry_run)
        
        cmd = 'python %s/index.py q.fst q.derep.fst --output q.index' %(self.library)
        self.ssub.submit_and_wait([cmd], self.dry_run)
        
        util.check_for_nonempty('q.index', dry_run=self.dry_run)
    
    def denovo_clustering(rename=True):
        '''Denovo clustering with USEARCH'''

        cmds = []
        for i in range(len(self.sids)):
            sid = self.sids[i]
            cmd = '%s -cluster_otus q.derep.fst -otus %s -otuid .%d' %(self.usearch, self.oi[i], sid)
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        util.check_for_nonempty(self.oi, self.dry_run)
        
        # Rename OTUs
        if rename == True:
            cmds = []
            for i in range(len(self.sids)):
                sid = self.sids[i]
                cmd = 'python %s/usearch_python/fasta_number.py %s OTU%d_ > %s' %(self.library, self.oi[i], sid, self.Oi[i])
                cmds.append(cmd)
            self.ssub.submit_and_wait(cmds, self.dry_run)
            util.check_for_nonempty(self.Oi, self.dry_run)
            self.ssub.move_files(self.Oi, self.oi, self.dry_run)
            util.check_for_nonempty(self.oi, self.dry_run)
    
    def remove_chimeras(self):
        '''Remove chimeras using gold database'''

        cmds = []
        for i in range(len(self.sids)):
            sid = self.sids[i]
            cmd = '%s -uchime_ref %s -db %s -nonchimeras %s -strand plus' %(self.usearch, self.oi[i], self.gold_db, self.Oi[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        
        util.check_for_nonempty(self.Oi, self.dry_run)
        self.ssub.move_files(self.Oi, self.oi, self.dry_run)
        util.check_for_nonempty(self.oi, self.dry_run)
    
    def reference_mapping(self):
        '''Map reads to reference databases'''
        
        util.check_for_nonempty(self.db, self.dry_run)
        util.check_for_collisions(self.uc, self.dry_run)

        cmds = []
        for i in range(len(self.sids)):
            cmd = '%s -usearch_global q.derep.fst -db %s -uc %s -strand both -id .%d' %(self.usearch, self.db[i], self.uc[i], self.reference_map_sids[i])
            cmds.append(cmd)
            
        self.ssub.submit_and_wait(cmds, self.dry_run)
        util.check_for_nonempty(self.uc, self.dry_run)
        
    def open_reference_mapping(self):
        '''Makes new otus from reads missed in the original reference mapping, then maps to those otus'''
        
        util.check_for_nonempty(self.uc, self.dry_run)
        util.check_for_collisions(self.open_fst, self.dry_run)        # unmatched seqs go into a new source fasta
        util.check_for_collisions(self.open_otu_tmp, self.dry_run)    # rep seqs for otus made from those unmatched seqs
        util.check_for_collisions(self.open_otu, self.dry_run)        # same rep seqs, but renamed as OTUs
        util.check_for_collisions(self.open_uc, self.dry_run)         # mapping information of unmatched seqs to denovo OTUs
        util.check_for_collisions(self.orig_uc, self.dry_run)         # mapping information for reference-based matching
        
        # grab the unmatched sequences from the database-reference uc
        cmds = ['python %s/uc2denovo.py q.derep.fst %s --output %s' %(self.library, self.uc[i], self.open_fst[i]) \
                for i in range(len(self.sids))]
        self.ssub.submit_and_wait(cmds, self.dry_run)
        util.check_for_nonempty(self.open_fst, self.dry_run)
        
        # denovo cluster those sequences
        cmds = ['%s -cluster_otus %s -otus %s -otu_radius_pct .%d' %(self.usearch, self.open_fst[i], self.open_otu[i], self.reference_map_pcts[i]) \
                for i in range(len(self.sids))]
        self.ssub.submit_and_wait(cmds, self.dry_run)
        util.check_for_nonempty(self.open_otu, self.dry_run)
        
        # rename the OTUs in those representative sequences
        cmds = ['python %s/usearch_python/fasta_number.py %s OTU%d_ > %s' %(self.library, self.open_otu[i], self.sids[i], self.open_otu_tmp[i]) \
                for i in range(len(self.sids))]
        self.ssub.submit_and_wait(cmds, self.dry_run)
        util.check_for_nonempty(self.open_otu_tmp, self.dry_run)
        
        # move the file with the renamed sequences into the spot as the final file
        self.ssub.move_files(self.open_otu_tmp, self.open_otu, self.dry_run)
        util.check_for_nonempty(self.open_otu, self.dry_run)
        
        # reference map against these new rep seqs
        cmds = ['%s -usearch_global %s -db %s -uc %s -strand both -id .%d' %(self.usearch, self.open_fst[i], self.open_otu[i], self.open_uc[i], self.reference_map_sids[i]) \
                for i in range(len(self.sids))]
        self.ssub.submit_and_wait(cmds, self.dry_run)
        util.check_for_nonempty(self.open_uc, self.dry_run)
        
        # move the original uc files to different place
        self.ssub.move_files(self.uc, self.orig_uc, self.dry_run)
        util.check_for_nonempty(self.orig_uc, self.dry_run)
        
        # combine the uc files
        cmds = ['grep "^H" %s > %s; cat %s >> %s' %(self.orig_uc[i], self.uc[i], self.open_uc[i], self.uc[i]) \
                for i in range(len(self.sids))]
        self.ssub.submit_and_wait(cmds, self.dry_run)
        util.check_for_nonempty(self.uc, self.dry_run)
    
    def make_otu_tables(self):
        '''Make OTU tables from UC file'''
        
        util.check_for_nonempty(self.uc + ['q.index'], self.dry_run)
        util.check_for_collisions(self.xi, self.dry_run)

        cmds = []
        for i in range(len(self.sids)):
            cmd = 'python %s/uc2otus.py %s q.index --output %s' %(self.library, self.uc[i], self.xi[i])
            
            # if we have a barcode file, use that order for the sample columns
            if self.barcodes is not None:
                cmd += ' --samples %s' %(self.barcodes)
            
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        
        util.check_for_nonempty(self.xi, self.dry_run)
        
    def make_seq_table(self):
        '''Make sequence table from the index file'''
        
        util.check_for_nonempty('q.fst', self.dry_run)
        util.check_for_collisions('seq.counts', self.dry_run)

        cmd = 'python %s/seq_table.py %s --output %s' %(self.library, 'q.fst', 'seq.counts')
        
        if self.barcodes is not None:
            cmd += ' --samples %s' %(self.barcodes)
        
        self.ssub.submit_and_wait([cmd], self.dry_run)
        util.check_for_nonempty('seq.counts', self.dry_run)
        
    def get_seq_tax(self):
        '''Get taxonomies for sequences in the seq table'''
        
        util.check_for_nonempty('seq.counts', self.dry_run)
        util.check_for_collisions(self.seq_tax_fn, self.dry_run)
        
        cmd = 'python %s/assign_seq_table_taxonomies.py %s --output %s' %(self.library, 'seq.counts', self.seq_tax_fn)
        
        self.ssub.submit_and_wait([cmd], self.dry_run)
        util.check_for_nonempty(self.seq_tax_fn, self.dry_run)
    

if __name__ == '__main__':
    # Initialize OTU caller
    oc = OTU_Caller()
    
    # Check fastq format
    if oc.check:
        message('Checking input formats')
        oc.check_format()
    
    # Split fastq
    if oc.split:
        message('Splitting fastq')
        oc.split_fastq()
        
    if oc.convert:
        message('Converting format')
        oc.convert_format()
    
    # Remove primers
    if oc.primers == True:
        message('Removing primers')
        oc.remove_primers()
    
    # Merge reads
    if oc.merge:
        message('Merging reads')
        oc.merge_reads()
        
    # Set current reads
    if oc.merge or oc.merged:
        oc.ci = oc.mi
    else:
        oc.ci = oc.fi
    
    # Demultiplex
    if oc.demultiplex == True:
        message('Demultiplexing')
        oc.demultiplex_reads()
    
    # Quality filter
    if oc.qfilter == True:
        message('Quality filtering')
        oc.quality_filter()
    
    # Dereplicate reads
    if oc.dereplicate == True:
        message('Dereplicating sequences')
        oc.dereplicate_reads()
        
    # Make index file
    if oc.index == True:
        message('Indexing samples')
        oc.make_index()
    
    # Denovo clustering
    if oc.denovo == True:
        message('Denovo clustering')
        oc.denovo_clustering(rename = True)
    
    # Map to reference database
    if oc.ref_gg == True:
        message('Mapping to reference')
        oc.reference_mapping()
        
    # Open reference clustering
    if oc.open_ref_gg:
        message('Open reference-based clustering')
        oc.open_reference_mapping()
        
    # Make sequence tables
    if oc.seq_table:
        message('Making sequence table')
        oc.make_seq_table()
        
    if oc.seq_tax:
        message('Assigning sequence table taxonomies')
        oc.get_seq_tax()
    
    # Make OTU tables
    if oc.otu_table == True:
        message('Making OTU tables')
        oc.make_otu_tables()
