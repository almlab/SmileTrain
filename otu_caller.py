#!/usr/bin/env python

'''
Main script in the pipeline. Produces lists of commands and submits them using ssub.
Options allow the user to run individual parts of the pipeline or the entire thing.
'''

import argparse, os, ConfigParser, subprocess, pickle
import ssub
import util
from util import *
import check_fastq_format

commands_fn = '.SmileTrain.commands.pkl'

class Submitter():
    '''runs jobs on a cluster, locally, or in a dry run'''
    def __init__(self, method, n_cpus=1, cluster=None):
        if method not in ['submit', 'local', 'dry_run']:
            raise ArgumentError("unexpcted method %s passed to Submitter" %(method))

        self.method = method

        if method == 'dry_run':
            self.dry_run = True
        elif method == 'submit':
            if cluster is None:
                raise RuntimeError("submitting jobs to cluster, but no cluster specified")

            self.dry_run = False

            # set up the cluster submission object
            self.ssub = ssub.Ssub(username=config.get('User', 'username'), cluster=config.get('User', 'cluster'),
                queue=config.get('User', 'queue'), tmp_dir=config.get('User', 'tmp_directory'), bashrc=config.get('Scripts', 'bashrc'), n_cpus=n_cpus)
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
            message("dry run: test that destinations are free: " + " ".join(fns), indent=4)
        else:
            util.check_for_collisions(fns)

    def check_is_executable(self, fn):
        '''check if a filename exists and is executable'''
        if not (os.path.isfile(fn) and os.access(fn, os.X_OK)):
            raise RuntimeError("file %s should be executable, but it is not" %(fn))

    def move_files(self, start_fns, end_fns):
        assert(len(start_fns) == len(end_fns))
        cmds = [['mv', x, y] for x, y in zip(start_fns, end_fns)]
        self.execute(cmds)

    def rm_files(self, fns):
        cmds = [['rm', fn] for fn in fns]
        self.execute(cmds)

    def execute(self, cmds):
        # recast all parts of the command to strings
        # swo> I regret this hack. It's to make spp's redirect command work.
        #cmds = [[str(x) for x in cmd] for cmd in cmds]
        def recast_cmd(cmd):
            if type(cmd) is str:
                return [cmd]
            elif type(cmd) is list:
                return [str(x) for x in cmd]

        cmds = [recast_cmd(cmd) for cmd in cmds]

        if self.method == 'submit':
            # recast commands as single lines
            cmds = [" ".join(cmd) for cmd in cmds]
            self.ssub.submit_and_wait(cmds)
        elif self.method == 'local':
            for cmd in cmds:
                print " ".join(cmd)
                if len(cmd) > 1:
                    subprocess.check_call(cmd)
                elif len(cmd) == 1:
                    subprocess.check_call(cmd[0], shell=True)
        elif self.method == 'dry_run':
            print "\n".join([" ".join(cmd) for cmd in cmds])


# open the config file sister to this script
config = ConfigParser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'user.cfg'))

def parse_args():
    '''arguments to be parsed and passed to the OTU_Caller object'''

    # create argument parser
    parser = argparse.ArgumentParser(description="Smile Train: a 16S pipeline", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # add groups
    group1 = parser.add_argument_group('Pipeline')
    group2 = parser.add_argument_group('Input files')
    group3 = parser.add_argument_group('Convert format')
    group4 = parser.add_argument_group('Remove primers')
    group5 = parser.add_argument_group('Merge reads')
    group6 = parser.add_argument_group('Demultiplexing Options')
    group7 = parser.add_argument_group('Quality filtering')
    group8 = parser.add_argument_group('Dereplicate')
    group9 = parser.add_argument_group('Chimeras')
    group10 = parser.add_argument_group('Indexing')
    group11 = parser.add_argument_group('Clustering')
    group12 = parser.add_argument_group('dbOTU options')
    group13 = parser.add_argument_group('Options')
    group14 = parser.add_mutually_exclusive_group()
    group_run = parser.add_mutually_exclusive_group()
    
    # add arguments
    group1.add_argument('--redo', action='store_true', help='Re-run with previously used options?')
    group1.add_argument('--all', action='store_true', help='Run primer, merge, demultiplex, filter, derep, index, ref_gg, and otu?')
    group1.add_argument('--check', action='store_true', help='Check input file format and intersection?')
    group1.add_argument('--split', action='store_true', help='Split the fastq files?')
    group1.add_argument('--convert', action='store_true', help='Convert fastq format?')
    group1.add_argument('--merge', action='store_true', help='Merge forward and reverse reads?')
    group1.add_argument('--primers', action='store_true', help='Remove primers?')
    group1.add_argument('--qfilter', default = False, action = 'store_true', help = 'Quality filter?')
    group1.add_argument('--ref_chimeras', action='store_true', help='Slay chimeras using reference database?')
    group1.add_argument('--chimeras', action='store_true', help='Slay chimeras de novo with usearch?')
    group1.add_argument('--dereplicate', action='store_true', help='Dereplicate?')
    group1.add_argument('--index', action='store_true', help='Make index file?')
    group1.add_argument('--denovo', default = False, action = 'store_true', help = 'Denovo clustering (UPARSE)?')
    group1.add_argument('--ref_gg', default = False, action = 'store_true', help = 'Reference mapping (Greengenes)?')
    group1.add_argument('--open_ref_gg', action='store_true', help='Reference map (Greengenes) and then denovo cluster?')
    group1.add_argument('--seq_table', action='store_true', help='Make a sequence table?')
    group1.add_argument('--seq_tax', action='store_true', help='Get taxonomies for sequence table?')
    group1.add_argument('--dbotu', action='store_true', help='Call OTUs using dbOTUs?')
    group1.add_argument('--otu_table', action='store_true', help='Make OTU table?')
    group2.add_argument('--forward', '-f', help='Input fastq (forward)')
    group2.add_argument('--reverse', '-r', help='Input fastq (reverse)')
    group2.add_argument('-p', help='Primer sequence (forward)')
    group2.add_argument('-q', help='Primer sequence (reverse)')
    group2.add_argument('--barcodes', '-b', default=None, help='Barcodes list')
    group4.add_argument('--p_mismatch', default=1, type=int, help='Number of mismatches allowed in primers')
    group6.add_argument('--b_mismatch', default=1, type=int, help='Number of mismatches allowed in barcodes')
    group7.add_argument('--truncqual', default = 2, type = int, help = '')
    group7.add_argument('--maxee', default = 2., type = float, help = 'Maximum expected error (UPARSE)')
    group7.add_argument('--trunclen', default=0, type=int, help='truncate all sequences to some length?')   # 0 means no truncation
    group9.add_argument('--gold_db', default=config.get('Data', 'gold'), help='Gold 16S database')
    group11.add_argument('--sids', default='91,94,97,99', help='Sequence identities for clustering')
    group12.add_argument('--alignref', default=config.get('dbOTU', 'alignref'), help='Reference alignment')
    group12.add_argument('--minlength', default=250, type=int, help='Minimum sequence length after alignment')
    group12.add_argument('--k_fold', default=0.0, type=float, help='k_fold change of OTU rep abundance over sequence to be merged')
    group12.add_argument('--pval', default=1e-4, type=float, help='p-value threshold for merge into existing OTU')
    group12.add_argument('--dbotu_chimeras', action='store_true', help='Remove chimeras de novo from dbOTUs?')
    group12.add_argument('--dbotu_split', action='store_true', help='Search for dbOTUs within pre-clustered sequences?')
    group12.add_argument('--dbotu_js', action='store_true', help='Merge statitically significantly different sequences if Jensen-Shannon divergence is below cut-off?')
    group12.add_argument('--dbotu_jscutoff', default=0.02, type=float, help='Jensen-Shannon divergence cut-off value used with --dbotu_js (default=0.02)')
    group12.add_argument('--dbotu_id', default=0.1, type=float, help='Distance used for dbOTUs and/or pre-clustering (default=0.1)')
    group13.add_argument('--align_start',  default=5, type=int, help='split upstream fastq into how many files?')
    group13.add_argument('--n_split', '-n', default=1, type=int, help='split upstream fastq into how many files?')
    group14.add_argument('--demultiplex', default = False, action = 'store_true', help = 'Demultiplex?')
    group14.add_argument('--already_demultiplexed', default = False, action = 'store_true', help = 'Already have seperate demultiplexed files?')
    group_run.add_argument('--dry_run', '-z', action='store_true', help='submit no jobs; suppress file checks; just print output commands')
    group_run.add_argument('--local', '-l', action='store_true', help='execute all tasks locally')
    
    # parse arguments
    if __name__ == '__main__':
        args = parser.parse_args()
    else:
        args = parser.parse_args('')

    if args.redo:
        # load command line arguments if available
        with open(commands_fn, 'rb') as f:
            args = pickle.load(f)
    else:
        # process arguments
        if args.all == True:
            args.check = args.split = args.convert = args.primers = args.merge = args.demultiplex = args.qfilter = args.dereplicate = args.index = args.ref_gg = args.otu_table = True
        args.sids = map(int, args.sids.split(','))
        
        # check combinations
        if args.check or args.split or args.convert or args.primers or args.merge or args.demultiplex or args.qfilter:
            if args.forward is None and args.reverse is None:
                raise RuntimeError("no fastq files selected")

        # save arguments for use with redo
        with open(commands_fn, 'wb') as f:
            pickle.dump(args, f)
        
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

        cluster = config.get('User', 'cluster')
        self.sub = Submitter(method, cluster=cluster, n_cpus=self.n_split)
    
    def get_filenames(self):
        '''Generate filenames to use in pipeline'''

        if self.forward or self.reverse:
            if self.forward:
                f_base = os.path.basename(self.forward)
            if self.reverse:
                r_base = os.path.basename(self.reverse)
                
            self.fi = ['%s.%d' %(self.forward, i) for i in range(self.n_split)] # forward reads (split)
            self.ri = ['%s.%d' %(self.reverse, i) for i in range(self.n_split)] # reverse reads (split)
            self.mi = ['%s.%d.merge' %(self.forward, i) for i in range(self.n_split)] # merged reads (split)
            self.Fi = ['%s.%d.tmp' %(self.forward, i) for i in range(self.n_split)] # forward reads (temp)
            self.Ri = ['%s.%d.tmp' %(self.reverse, i) for i in range(self.n_split)] # reverse reads (temp)
            self.Mi = ['%s.%d.tmp' %(self.forward, i) for i in range(self.n_split)] # merged reads (temp)
            self.ci = ['q.%d.fst' %(i) for i in range(self.n_split)] # current reads
            self.Ci = ['q.%d.tmp' %(i) for i in range(self.n_split)] # current reads (temp)
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

        if self.forward and self.reverse:
            cmds += ['python', '%s/check_intersect.py', self.forward, self.reverse]

        self.sub.execute(cmds)
    
    def split_fastq(self):
        '''Split forward and reverse reads (for parallel processing)'''

        # do forward only if there is a forward file; similar for reverse
        do_forward = self.forward
        do_reverse = self.reverse

        # check for inputs and collisions of output
        if do_forward:
            self.sub.check_for_nonempty(self.forward)
            self.sub.check_for_collisions(['%s.%s' %(self.forward, i) for i in range(self.n_split)])
        if do_reverse:
            self.sub.check_for_nonempty(self.reverse)
            self.sub.check_for_collisions(['%s.%s' %(self.reverse, i) for i in range(self.n_split)])
        
        # Get list of commands
        cmds = []
        if do_forward:
            cmd = ['python', '%s/split_fastq.py' %(self.library), self.forward, self.n_split]
            cmds.append(cmd)
        if do_reverse:
            cmd = ['python', '%s/split_fastq.py' %(self.library), self.reverse, self.n_split]
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
        for i in range(self.n_split):
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

    def merge_reads(self):
        '''Merge forward and reverse reads using USEARCH'''

        # check for inputs and collisions
        self.sub.check_for_nonempty(self.fi + self.ri)
        self.sub.check_for_collisions(self.Fi + self.Ri)
        
        # check that usearch is ready to go
        self.sub.check_is_executable(self.usearch)
        
        # check that forward and reverse reads intersect
        cmds = []
        for i in range(self.n_split):
            cmd = ['python', '%s/check_intersect.py' %(self.library), self.fi[i], self.ri[i]]
            cmds.append(cmd)
        self.sub.execute(cmds)
        
        # Merge reads
        cmds = []
        for i in range(self.n_split):
            cmd = [self.usearch, '-fastq_mergepairs', self.fi[i], '-reverse', self.ri[i], '-fastq_truncqual', self.truncqual, '-fastqout', self.Fi[i]]
            cmds.append(cmd)
        self.sub.execute(cmds)
        
        self.sub.check_for_nonempty(self.Fi)
        self.sub.rm_files(self.fi + self.ri)
        self.sub.move_files(self.Fi, self.fi)
        self.sub.check_for_nonempty(self.fi)
    
    def remove_primers(self):
        '''Remove diversity region + primer and discard reads with > 2 mismatches'''

        # do forward only if there is a forward read file and a forward primer
        # similar for reverse
        both = self.p and self.q
        if not both:
            if self.p and not self.q:
                forward_only = True
            else:
                raise RuntimeError("remove primers called with bad input: need -p or both -p and -q")

        # check for inputs and collisions of output
        self.sub.check_for_nonempty(self.fi)
        self.sub.check_for_collisions(self.Fi)
        
        # get list of commands using forward only
        cmds = [['python', '%s/remove_primers.py' %(self.library), fi, self.p, '--max_primer_diffs', self.p_mismatch, '--output', fo] for fi, fo in zip(self.fi, self.Fi)]

        # add reverse primer if needed
        if both:
            cmds = [cmd + ['--reverse_primer', self.q] for cmd in cmds]
        
        # submit commands
        self.sub.execute(cmds)

        # validate output and move files
        self.sub.check_for_nonempty(self.Fi)
        self.sub.move_files(self.Fi, self.fi)
        self.sub.check_for_nonempty(self.fi)
    
    def demultiplex_reads(self):
        '''Demultiplex samples using index and barcodes'''
        
        self.sub.check_for_nonempty(self.ci)
        self.sub.check_for_collisions(self.Ci)

        cmds = []
        for i in range(self.n_split):
            cmd = ['python', '%s/map_barcodes.py' %(self.library), self.ci[i], self.barcodes, '--max_barcode_diffs', self.b_mismatch, '--output', self.Ci[i]]
            cmds.append(cmd)
        self.sub.execute(cmds)
        
        self.sub.check_for_nonempty(self.Ci)
        self.sub.move_files(self.Ci, self.ci)
        self.sub.check_for_nonempty(self.ci)

    def reformat_headers(self):
        self.sub.check_for_nonempty(self.ci)
        self.sub.check_for_collisions(self.Ci)

        cmds = []
        for i in range(self.n_split):
            cmd = ['python', '%s/reformat_headers.py' %(self.library), self.ci[i],'--output', self.Ci[i]]
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
        self.sub.check_is_executable(self.usearch)
        
        # check that the files are in the right format
        if self.dry_run:
            cmds = ['python %s/check_fastq_format.py %s' %(self.library, f) for f in self.ci]
            message("\n".join(cmds), indent=0)
        else:
            check_fastq_format.check_illumina_format(self.ci, ['illumina18', 'ambiguous'])

        cmds = []
        for i in range(self.n_split):
            cmd = [self.usearch, '-fastq_filter', self.ci[i], '-fastq_truncqual', self.truncqual, '-fastq_maxee', self.maxee, '-fastaout', self.Ci[i]]

            if self.trunclen > 0:
                cmd += ['-fastq_trunclen', self.trunclen]

            cmds.append(cmd)
        self.sub.execute(cmds)
        
        cmd = ['python', '%s/combine_fasta.py' %(self.library), '--output', 'q.fst'] + self.Ci
        #swo> I regret this hack; i need to [] the cmd
        self.sub.execute([cmd])
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
    
    def dbotu_progressive_clustering(self):
       '''Denovo clustering with USEARCH'''
       perllib = config.get('dbOTU', 'perllib')
       cmds = []
       cmds.append(['perl', '%s/find_replace_seq_dash-period.pl' % perllib, 'unique.good.align', 'unique.good.align.ng'])
       cmds.append(['perl', '%s/fasta2uchime_size.pl' % perllib, 'unique.f0.good.mat', 'unique.good.align.ng', 'unique.good.align.ng.size'])
       cmds.append([self.usearch, '-cluster_otus', 'unique.good.align.ng.size', '-uparseout', 'unique.97.uc', '-otus', 'unique.97.otus.fa'])
       cmds.append(['perl', '%s/UC2list4.pl' % perllib, 'unique.97.uc', 'unique.97.uc.list'])
       cmds.append([self.usearch, '-sortbylength', 'unique.97.otus.fa', '-fastaout', 'unique.97.sorted.fa'])
       uc_list=['unique.97.uc.list']
       upper=96
       lower=int(100*(1-self.dbotu_id)) - 1
       for i in range(upper, lower, -1):
           previous= i + 1
           cmds.append([self.usearch, '-cluster_smallmem', 'unique.%d.sorted.fa' % previous, '-id', '0.%d' % i, '--uc', 'unique.%d.uc' % i, '-centroids', 'unique.%d.otus.fa' % i])
           cmds.append(['perl', '%s/UC2list3.pl' % perllib, 'unique.%d.uc' % i, 'unique.%d.uc.list' % i])
           cmds.append([self.usearch, '-sortbylength', 'unique.%d.otus.fa' % i, '-fastaout', 'unique.%d.sorted.fa' % i])
           uc_list.append(str('unique.%d.uc.list' % i))
       input_lists=",".join([str(x) for x in  uc_list])
       cmds.append(['perl', '%s/merge_progressive_clustering4.pl'% perllib, input_lists, 'unique.PC.final.list'])

       self.sub.execute(cmds)
       self.sub.check_for_nonempty('unique.PC.final.list')       
 
    def remove_reference_chimeras(self):
        '''Remove chimeras using gold database'''
        cmds = []
        for i in range(len(self.sids)):
            sid = self.sids[i]
            cmd = [self.usearch, '-uchime_ref', self.oi[i], '-db', self.gold_db, '-nonchimeras', self.Oi[i], '-strand', 'plus']
            cmds.append(cmd)
        self.sub.execute(cmds)
        
        self.sub.check_for_nonempty(self.Oi)
        self.sub.move_files(self.Oi, self.oi)
        self.sub.check_for_nonempty(self.oi)

    def dbotu_remove_chimeras(self):
        '''remove dbotu chimeras de novo'''
        perllib = config.get('dbOTU', 'perllib')

        cmds = []
        cmd = 'perl %s/fasta2filter_from_mat_SmileTrain.pl unique.dbOTU.mat q.derep.fst > unique.dbOTU.ng.fasta' % perllib
        cmds.append(cmd)
        cmd = '%s -uchime_denovo unique.dbOTU.ng.fasta -nonchimeras unique.dbOTU.nonchimera.fasta' % self.usearch
        cmds.append(cmd)
        cmd = 'perl %s/filter_mat_from_fasta_SmileTrain.pl unique.dbOTU.mat unique.dbOTU.nonchimera.fasta > unique.dbOTU.nonchimera.mat' % perllib
        cmds.append(cmd)

        self.sub.execute(cmds)
        self.sub.check_for_nonempty('unique.dbOTU.nonchimera.fasta')
    
    def remove_denovo_chimeras(self):
        '''remove chimeras identified de novo by usearch'''
        cmds = [[self.usearch, '-uchime_denovo', self.oi[i], '-nonchimeras', self.Oi[i], '-strand', 'plus'] for i in range(len(self.sids))]
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

    def dbotu_alignment(self):
        '''Call otus using dbotu algorithm'''
        perllib = config.get('dbOTU', 'perllib')
        mothur = config.get('dbOTU', 'mothur')
        caller = config.get('dbOTU', 'caller')

        cmds = []
        cmds.append(['perl', '%s/temp_071514.pl' % perllib, 'q.derep.fst', 'q.index', 'unique'])
        cmds.append(['%s "#align.seqs(fasta=unique.fa, reference=%s)"' %(mothur, self.alignref)])
        cmds.append(['%s "#screen.seqs(fasta=unique.align, start=%d, minlength=%d)"' %(mothur, self.align_start, self.minlength)])
        cmds.append('perl %s/filter_mat_from_fasta.pl unique.f0.mat unique.good.align > unique.f0.good.mat' %(perllib))

        self.sub.execute(cmds)
        self.sub.check_for_nonempty(['unique.good.align', 'unique.f0.good.mat'])
        
    def dbotu_call_otus(self):
        '''Remove redundancy and errors with dbotus'''
        mothur = config.get('dbOTU', 'mothur')
        caller = config.get('dbOTU', 'caller')
        cmds = []
        dbcmd = ['python', caller, 'unique.f0.good.mat', 'unique.good.align', 'unique.dbOTU', '-k', str(self.k_fold), '-p', str(self.pval), '-d', str(self.dbotu_id)]
        if oc.dbotu_split:
            dbcmd.append('-s')
            dbcmd.append('unique.PC.final.list')
        if oc.dbotu_js:
            dbcmd.append('--useJS')
            dbcmd.append(str(self.dbotu_jscutoff))

        
        cmds.append(dbcmd)
        cmds.append(['%s "#degap.seqs(fasta=unique.dbOTU.fasta)"' %(mothur)])

        self.sub.execute(cmds)

        self.sub.check_for_nonempty(['unique.dbOTU.list', 'unique.dbOTU.ng.fasta', 'unique.dbOTU.mat', 'unique.dbOTU.log'])
    
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

    # Merge reads
    if oc.merge:
        message('Merging reads')
        oc.merge_reads()
    
    # Remove primers
    if oc.primers == True:
        message('Removing primers')
        oc.remove_primers()
        
    # Set current reads
    # swo> obsolete, now that there are no separate merged files
    if hasattr(oc, 'fi'):
        oc.ci = oc.fi
    
    # Demultiplex
    if oc.demultiplex == True:
        message('Demultiplexing')
        oc.demultiplex_reads()

    if oc.already_demultiplexed == True:
        message('Already demultiplexed: reformatting sequence headers')
        oc.reformat_headers()
    
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
    if oc.ref_gg or oc.open_ref_gg:
        message('Mapping to reference')
        oc.reference_mapping()

    # Call dbOTUs
    if oc.dbotu:
        message("dbOTU: aligning sequences")
        oc.dbotu_alignment()

        if oc.dbotu_split:
            message("dbOTU: progressive clustering")
            oc.dbotu_progressive_clustering()

        message("Calling dbOTUs")
        oc.dbotu_call_otus()

    # Chimera removal
    if oc.ref_chimeras:
        message("Removing chimeras by reference")
        oc.remove_reference_chimeras()

    if oc.chimeras:
        message("Removing chimeras de novo with uchime")
        oc.remove_denovo_chimeras()

    if oc.dbotu_chimeras:
        message("Removing chimeras from dbOTUs de novo")
        oc.dbotu_remove_chimeras()
        
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
