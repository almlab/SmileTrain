#!/usr/bin/env python

'''
Main script in the pipeline. Produces lists of commands and submits them using ssub.
Options allow the user to run individual parts of the pipeline or the entire thing.
'''

import argparse, os, ConfigParser, subprocess
import ssub
import util
from util import *

class Submitter():
    '''runs jobs on a cluster, locally, or in a dry run'''
    def __init__(self, method, n_cpus=1):
        if method not in ['submit', 'local', 'dry_run']:
            raise ArgumentError("unexpcted method %s passed to Submitter" %(method))

        self.method = method

        if method == 'dry_run':
            self.dry_run = True
        elif method == 'submit':
            self.dry_run = False
            self.ssub = ssub.Ssub()
            self.ssub.n_cpus = n_cpus
        elif method == 'local':
            self.dry_run = False

    def check_for_existence(self, fns):
        '''assert that each of filenames does exist'''

        fns = util.listify(fns)

        if self.method == 'dry_run':
            message("dry run: test for existence of files: " + " ".join(fns), indent=4)
        else:
            # running ls first seems to prevent spurious empties
            subprocess.check_output(['ls', '-lah'])
            tests = [os.path.isfile(fn) for fn in fns]
            if False in tests:
                bad_names = " ".join([fn for fn, test in zip(fns, tests) if test == False])
                raise RuntimeError("file(s) missing: %s" % bad_names)

    def check_for_nonempty(self, fns):
        '''assert that each file exists and is nonempty'''
        fns = util.listify(fns)
    
        if self.method == 'dry_run':
            message("dry run: test that files are non-empty: " + " ".join(fns), indent=4)
        else:
            self.check_for_existence(fns)
            tests = [os.stat(fn).st_size > 0 for fn in fns]
            if False in tests:
                bad_names = " ".join([fn for fn, t in zip(fns, tests) if t == False])
                raise RuntimeError("file(s) empty: " + bad_names)

    def check_for_collisions(self, fns):
        if self.method == 'dry_run':
            message("dry run: test that destinations are free: " + " ".join(filenames), indent=4)
        else:
            util.check_for_collisions(fns)

    def check_is_executable(self, fn):
        '''check if a filename exists and is executable'''
        if not (os.path.isfile(fn) and os.access(fn, os.X_OK)):
            raise RuntimeError("file %s should be executable, but it is not" %(fn))

    def move_files(self, start_fns, end_fns):
        cmds = [['mv', x, y] for x, y in zip(start_fns, end_fns)]
        self.execute(cmds)

    def rm_files(self, fns):
        cmds = [['rm', fn] for fn in fns]
        self.execute(cmds)

    def execute(self, cmds):
        # recast all parts of the command to strings
        cmds = [[str(x) for x in cmd] for cmd in cmds]

        if self.method == 'submit':
            self.ssub.submit_and_wait(cmds)
        elif self.method == 'local':
            for cmd in cmds:
                print " ".join(cmd)
                subprocess.check_call(cmd)
        elif self.method == 'dry_run':
            print "\n".join([" ".join(cmd) for cmd in cmds])


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
    group_run = parser.add_mutually_exclusive_group()
    
    # add arguments
    group1.add_argument('--all', action='store_true', help='Run primer, merge, demultiplex, filter, derep, index, ref_gg, and otu?')
    group1.add_argument('--check', action='store_true', help='Check input file format?')
    group1.add_argument('--split', action='store_true', help='Split the fastq files?')
    group1.add_argument('--convert', action='store_true', help='Convert fastq format?')
    group1.add_argument('--primers', action='store_true', help='Remove primers?')
    group1.add_argument('--merge', action='store_true', help='Merge forward and reverse reads?')
    group1.add_argument('--demultiplex', default = False, action = 'store_true', help = 'Demultiplex?')
    group1.add_argument('--qfilter', default = False, action = 'store_true', help = 'Quality filter?')
    group1.add_argument('--ref_chimera', action='store_true', help='Slay chimeras using reference database?')
    group1.add_argument('--chimera', action='store_true', help='Slay chimeras de novo with usearch?')
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
    group7.add_argument('--trunclen', default=0, type=int, help='truncate all sequences to some length?')   # 0 means no truncation
    group9.add_argument('--gold_db', default=config.get('Data', 'gold'), help='Gold 16S database')
    group11.add_argument('--sids', default='91,94,97,99', help='Sequence identities for clustering')
    group12.add_argument('--n_cpus', '-n', default = 1, type = int, help='Number of CPUs')
    group_run.add_argument('--dry_run', '-z', action='store_true', help='submit no jobs; suppress file checks; just print output commands')
    group_run.add_argument('--local', '-l', action='store_true', help='execute all tasks locally')
    
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
        if self.dry_run:
            method = 'dry_run'
        elif self.local:
            method = 'local'
        else:
            method = 'submit'

        self.sub = Submitter(method, n_cpus=self.n_cpus)
    
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
        
        self.open_fst = ['q.%d.no_match.fst' %(sid) for sid in self.sids] # unmatched seqs from ref mapping
        self.seq_tax_fn = 'seq.tax'
        
        # if reference-based clustering at 97%, make sure reads are 98.5% dissimilar
        self.reference_map_sids = [0.5*(100.0+float(sid)) for sid in self.sids]
        self.reference_map_pcts = [100.0 - sid for sid in self.reference_map_sids]
        
        # Get database for read mapping
        if self.denovo == True:
            self.db = self.oi
        elif self.ref_gg or self.open_ref_gg:
            self.db = ['%s/%d_otus.fasta' %(self.ggdb, sid) for sid in self.sids]
            
    def check_format(self):
        '''Make sure we have the correct input format'''
        files = []
        if self.forward:
            files.append(self.forward)
        
        if self.reverse:
            files.append(self.reverse)
            
        message('Testing format of %s' %(" ".join(files)))

        cmds = [['python', '%s/check_fastq_format.py' %(self.library), f] for f in files]
        self.sub.execute(cmds)
    
    def split_fastq(self):
        '''Split forward and reverse reads (for parallel processing)'''

        # do forward only if there is a forward file; similar for reverse
        do_forward = self.forward
        do_reverse = self.reverse

        # check for inputs and collisions of output
        if do_forward:
            self.sub.check_for_nonempty(self.forward)
            self.sub.check_for_collisions(['%s.%s' %(self.forward, i) for i in range(self.n_cpus)])
        if do_reverse:
            self.sub.check_for_nonempty(self.reverse)
            self.sub.check_for_collisions(['%s.%s' %(self.reverse, i) for i in range(self.n_cpus)])
        
        # Get list of commands
        cmds = []
        if do_forward:
            cmd = ['python', '%s/split_fastq.py' %(self.library), self.forward, self.n_cpus]
            cmds.append(cmd)
        if do_reverse:
            cmd = ['python', '%s/split_fastq.py' %(self.library), self.reverse, self.n_cpus]
            cmds.append(cmd)
        
        # submit commands
        self.sub.execute(cmds)

        # validate output
        if do_forward:
            self.sub.check_for_nonempty(self.fi)
        if do_reverse:
            self.sub.check_for_nonempty(self.ri)
            
    def convert_format(self):
        '''Convert to compatible fastq format'''
        
        if self.forward:
            self.sub.check_for_nonempty(self.fi)
            self.sub.check_for_collisions(self.Fi)
            
        if self.reverse:
            self.sub.check_for_nonempty(self.ri)
            self.sub.check_for_collisions(self.Fi)
            
        cmds = []
        for i in range(self.n_cpus):
            if self.forward:
                cmd = ['python', '%s/convert_fastq.py' %(self.library), self.fi[i], '--output', self.Fi[i]]
                cmds.append(cmd)
            if self.reverse:
                cmd = ['python', '%s/convert_fastq.py' %(self.library), self.ri[i], '--output', self.Ri[i]]
                cmds.append(cmd)
                
        self.sub.execute(cmds)
        
        # validate output and move files
        if self.forward:
            self.sub.check_for_nonempty(self.Fi)
            self.sub.move_files(self.Fi, self.fi)
            self.sub.check_for_nonempty(self.fi)

        if self.reverse:
            self.sub.check_for_nonempty(self.Ri)
            self.sub.move_files(self.Ri, self.ri)
            self.sub.check_for_nonempty(self.ri)
    
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
            self.sub.check_for_nonempty(self.fi)
            self.sub.check_for_collisions(self.Fi)
        if do_reverse:
            self.sub.check_for_nonempty(self.ri)
            self.sub.check_for_collisions(self.Ri)
        
        # get list of commands using forward, reverse, or both
        cmds = []
        for i in range(self.n_cpus):
            if do_forward:
                cmd = ['python', '%s/remove_primers.py' %(self.library), self.fi[i], self.p, '--max_primer_diffs', self.p_mismatch, '--output', self.Fi[i]]
                cmds.append(cmd)
            if do_reverse:
                cmd = ['python', '%s/remove_primers.py' %(self.library), self.ri[i], self.p, '--max_primer_diffs', self.p_mismatch, '--output', self.Ri[i]]
                cmds.append(cmd)
        
        # submit commands
        self.sub.execute(cmds)

        # validate output and move files
        if do_forward:
            self.sub.check_for_nonempty(self.Fi)
            self.sub.move_files(self.Fi, self.fi)
            self.sub.check_for_nonempty(self.fi)

        if do_reverse:
            self.sub.check_for_nonempty(self.Ri)
            self.sub.move_files(self.Ri, self.ri)
            self.sub.check_for_nonempty(self.ri)
    
    def merge_reads(self):
        '''Merge forward and reverse reads using USEARCH'''

        # check for inputs and collisions
        self.sub.check_for_nonempty(self.fi + self.ri)
        self.sub.check_for_collisions(self.Fi + self.Ri)
        
        # check that usearch is ready to go
        self.sub.check_is_executable(self.usearch)
        
        # check that forward and reverse reads intersect
        cmds = []
        for i in range(self.n_cpus):
            cmd = ['python', '%s/check_intersect.py' %(self.library), self.fi[i], self.ri[i]]
            cmds.append(cmd)
        self.sub.execute(cmds)
        
        # Merge reads
        cmds = []
        for i in range(self.n_cpus):
            cmd = [self.usearch, '-fastq_mergepairs', self.fi[i], '-reverse', self.ri[i], '-fastq_truncqual', self.truncqual, '-fastqout', self.mi[i]]
            cmds.append(cmd)
        self.sub.execute(cmds)
        
        self.sub.check_for_nonempty(self.mi)
        self.sub.rm_files(self.fi + self.ri)
    
    def demultiplex_reads(self):
        '''Demultiplex samples using index and barcodes'''
        
        self.sub.check_for_nonempty(self.ci)
        self.sub.check_for_collisions(self.Ci)

        cmds = []
        for i in range(self.n_cpus):
            cmd = ['python', '%s/map_barcodes.py' %(self.library), self.ci[i], self.barcodes, '--max_barcode_diffs', self.b_mismatch, '--output', self.Ci[i]]
            cmds.append(cmd)
        self.sub.execute(cmds)
        
        self.sub.check_for_nonempty(self.Ci)
        self.sub.move_files(self.Ci, self.ci)
        self.sub.check_for_nonempty(self.ci)
    
    def quality_filter(self):
        '''Quality filter with truncqual and maximum expected error'''
        
        # validate input/output
        self.sub.check_for_nonempty(self.ci)
        self.sub.check_for_collisions(self.Ci)
        
        # check that usearch is ready to go
        assert(self.sub.is_executable(self.usearch))
        
        # check that the files are in the right format
        if self.dry_run:
            cmds = ['python %s/check_fastq_format.py %s' %(self.library, f) for f in self.ci]
            message("\n".join(cmds), indent=0)
        else:
            check_fastq_format.check_illumina_format(self.ci, ['illumina18', 'ambiguous'])

        cmds = []
        for i in range(self.n_cpus):
            cmd = [self.usearch, '-fastq_filter', self.ci[i], '-fastq_truncqual', self.truncqual, '-fastq_maxee', self.maxee, '-fastaout', self.Ci[i]]

            if args.trunclen > 0:
                cmd += ['-fastq_trunclen', str(args.trunclen)]

            cmds.append(cmd)
        self.sub.submit_and_wait(cmds)
        
        cmd = ['python', '%s/combine_fasta.py' %(self.library), '--output', 'q.fst'] + self.Ci
        self.sub.execute(cmd)
        self.sub.check_for_nonempty('q.fst')
        self.sub.rm_files(self.Ci)
    
    def dereplicate_reads(self):
        '''Concatenate files and dereplicate'''

        cmd = ['python',  '%s/derep_fulllength.py' %(self.library), 'q.fst', '--output', 'q.derep.fst']
        self.sub.execute([cmd])
        self.sub.check_for_nonempty('q.derep.fst')
        
    def make_index(self):
        '''Make an index file'''
        # verify input & check for collisions with output
        self.sub.check_for_nonempty(['q.fst', 'q.derep.fst'])
        self.sub.check_for_collisions('q.index')
        
        cmd = ['python', '%s/index.py' %(self.library), 'q.fst', 'q.derep.fst', '--output', 'q.index']
        self.sub.execute([cmd])
        
        self.sub.check_for_nonempty('q.index')
    
    def denovo_clustering(rename=True):
        '''Denovo clustering with USEARCH'''
        cmds = []
        for i in range(len(self.sids)):
            sid = self.sids[i]
            cmd = [self.usearch, '-cluster_otus', 'q.derep.fst', '-otus', self.oi[i], '-otuid', '.%d' %(sid)]
            cmds.append(cmd)
        self.sub.execute(cmds)
        self.sub.check_for_nonempty(self.oi)
        
        # Rename OTUs
        if rename == True:
            cmds = []
            for i in range(len(self.sids)):
                sid = self.sids[i]
                cmd = ['python', '%s/usearch_python/fasta_number.py' %(self.library), self.oi[i], 'OTU%d_' %(sid), '>', self.Oi[i]]
                cmds.append(cmd)
            self.sub.execute(cmds)
            self.sub.check_for_nonempty(self.Oi)
            self.sub.move_files(self.Oi, self.oi)
            self.sub.check_for_nonempty(self.oi)
    
    def remove_reference_chimeras(self):
        '''Remove chimeras using gold database'''
        cmds = []
        for i in range(len(self.sids)):
            sid = self.sids[i]
            cmd = [self.usearch, '-uchime_ref', self.oi[i], '-db', self.gold_db, '-nochimeras', self.Oi[i], '-strand', 'plus']
            cmds.append(cmd)
        self.sub.execute(cmds)
        
        self.sub.check_for_nonempty(self.Oi)
        self.sub.move_files(self.Oi, self.oi)
        self.sub.check_for_nonempty(self.oi)
    
    def remove_denovo_chimeras(self):
        '''remove chimeras identified de novo by usearch'''
        cmds = [[self.usearch, '-uchime_denovo', self.oi[i], '-nochimeras', self.Oi[i], '-strand', 'plus'] for i in range(len(self.sids))]
        self.sub.execute(cmds)

        self.sub.check_for_nonempty(self.Oi)
        self.sub.move_files(self.Oi, self.oi)
        self.sub.check_for_nonempty(self.oi)

    def reference_mapping(self):
        '''Map reads to reference databases'''
        self.sub.check_for_nonempty(self.db)
        self.sub.check_for_collisions(self.uc)

        cmds = []
        for i in range(len(self.sids)):
            cmd = [self.usearch, '-usearch_global', 'q.derep.fst', '-db', self.db[i], '-uc', self.uc[i], '-strand', 'both', '-id', '.%d' %(self.reference_map_sids[i])]

            if self.open_ref_gg:
                cmd += '-notmatched', self.open_fst[i]

            cmds.append(cmd)
            
        self.sub.execute(cmds)
        self.sub.check_for_nonempty(self.uc)

        if self.open_ref_gg:
            message("unmatches sequences left in following files. move them to a new work folder for de novo clustering.")
            print "\n".join(self.open_fst)
    
    def make_otu_tables(self):
        '''Make OTU tables from uc file'''    
        self.sub.check_for_nonempty(self.uc + ['q.index'])
        self.sub.check_for_collisions(self.xi)

        cmds = []
        for i in range(len(self.sids)):
            cmd = ['python', '%s/uc2otus.py' %(self.library), self.uc[i], 'q.index', '--output', self.xi[i]]
            
            # if we have a barcode file, use that order for the sample columns
            if self.barcodes is not None:
                cmd += ['--samples', self.barcodes]
            
            cmds.append(cmd)
        self.sub.execute(cmds)
        
        self.sub.check_for_nonempty(self.xi)
        
    def make_seq_table(self):
        '''Make sequence table from the index file'''
        self.sub.check_for_nonempty('q.fst')
        self.sub.check_for_collisions('seq.counts')

        cmd = ['python', '%s/seq_table.py' %(self.library), 'q.fst', 'q.derep.fst', '--output', 'seq.counts']
        
        if self.barcodes is not None:
            cmd += ['--samples', self.barcodes]
        
        self.sub.execute([cmd])
        self.sub.check_for_nonempty('seq.counts')
        
    def get_seq_tax(self):
        '''Get taxonomies for sequences in the seq table'''
        self.sub.check_for_nonempty('seq.counts')
        self.sub.check_for_collisions(self.seq_tax_fn)
        
        cmd = ['python', '%s/assign_seq_table_taxonomies.py' %(self.library), 'seq.counts', '--output', self.seq_tax_fn]
        
        self.sub.execute([cmd])
        self.sub.check_for_nonempty(self.seq_tax_fn)
    

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

    # Chimera removal
    if oc.ref_chimera:
        message("Removing chimeras by reference")
        oc.remove_reference_chimeras()

    if oc.chimera:
        message("Removing chimeras de novo with uchime")
        oc.remove_denovo_chimeras()
    
    # Map to reference database
    if oc.ref_gg or oc.open_ref_gg:
        message('Mapping to reference')
        oc.reference_mapping()
        
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