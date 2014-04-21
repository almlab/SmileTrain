#!/usr/bin/env python

'''
Main script in the pipeline. Produces lists of commands and submits them using ssub.
Options allow the user to run individual parts of the pipeline or the entire thing.
'''

import argparse, os, ConfigParser
import ssub
from util import *

# open the config file sister to this script
config = ConfigParser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'train.cfg'))

def parse_args():
    '''arguments to be parsed and passed to the OTU_Caller object'''

    # create argument parser
    parser = argparse.ArgumentParser()
    
    # add groups
    group1 = parser.add_argument_group('Pipeline')
    group2 = parser.add_argument_group('Input files')
    group3 = parser.add_argument_group('Remove primers')
    group4 = parser.add_argument_group('Merge reads')
    group5 = parser.add_argument_group('Demultiplex')
    group6 = parser.add_argument_group('Quality filtering')
    group7 = parser.add_argument_group('Dereplicate')
    group8 = parser.add_argument_group('Chimeras')
    group9 = parser.add_argument_group('Clustering')
    group10 = parser.add_argument_group('Options')
    
    # add arguments
    group1.add_argument('--all', default = False, action = 'store_true', help = 'Run all steps of pipeline?')
    group1.add_argument('--dont_split', action='store_true', help='Fastq files already split?')
    group1.add_argument('--primers', default = False, action = 'store_true', help = 'Remove primers?')
    group1.add_argument('--merge', default = False, action = 'store_true', help = 'Merge forward and reverse reads?')
    group1.add_argument('--demultiplex', default = False, action = 'store_true', help = 'Demultiplex?')
    group1.add_argument('--qfilter', default = False, action = 'store_true', help = 'Quality filter?')
    group1.add_argument('--chimeras', default = False, action = 'store_true', help = 'Chimera slay?')
    group1.add_argument('--dereplicate', action='store_true', help='Dereplicate?')
    group1.add_argument('--denovo', default = False, action = 'store_true', help = 'Denovo clustering (UPARSE)?')
    group1.add_argument('--ref_gg', default = False, action = 'store_true', help = 'Reference mapping (Greengenes)?')
    group1.add_argument('--otu_table', action='store_true', help='Make OTU table?')
    group2.add_argument('-f', help = 'Input fastq (forward)')
    group2.add_argument('-r', help = 'Input fastq (reverse)')
    group2.add_argument('-p', help = 'Primer sequence (forward)')
    group2.add_argument('-q', help = 'Primer sequence (reverse)')
    group2.add_argument('-b', help = 'Barcodes list')
    group2.add_argument('-x', help = 'Index fastq')
    group3.add_argument('--p_mismatch', default=1, type=int, help='Number of mismatches allowed in primers')
    group5.add_argument('--b_mismatch', default=1, type=int, help='Number of mismatches allowed in barcodes')
    group6.add_argument('--truncqual', default = 2, type = int, help = '')
    group6.add_argument('--maxee', default = 2., type = float, help = 'Maximum expected error (UPARSE)')
    group8.add_argument('--gold_db', default=config.get('Data', 'gold'), help='Gold 16S database')
    group9.add_argument('--sids', default='91,94,97,99', help='Sequence identities for clustering')
    group10.add_argument('--n_cpus', '-n', default = 1, type = int, help='Number of CPUs')
    group10.add_argument('--dry_run', '-z', action='store_true', help='Submit no jobs; just print output commands')
    
    # parse arguments
    if __name__ == '__main__':
        args = parser.parse_args()
    else:
        args = parser.parse_args('')
    
    # process arguments
    if args.all == True:
        args.primers = args.merge = args.demultiplex = args.qfilter = args.chimeras = args.ref_gg = True
    args.sids = map(int, args.sids.split(','))
        
    return args


class OTU_Caller():
    '''
    A namespace and method container for the parsed command line options. Submitted jobs
    refer to other scripts in the library.
    '''

    def __init__(self):
        # initialize variables
        self.usearch = 'usearch'
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

        if self.f:
            f_base = os.path.basename(self.f)
        if self.r:
            r_base = os.path.basename(self.r)
        self.fi = ['%s.%d' %(self.f, i) for i in range(self.n_cpus)] # forward reads (split)
        self.ri = ['%s.%d' %(self.r, i) for i in range(self.n_cpus)] # reverse reads (split)
        self.mi = ['%s.%d.merge' %(self.f, i) for i in range(self.n_cpus)] # merged reads (split)
        self.Fi = ['%s.%d.tmp' %(self.f, i) for i in range(self.n_cpus)] # forward reads (temp)
        self.Ri = ['%s.%d.tmp' %(self.r, i) for i in range(self.n_cpus)] # reverse reads (temp)
        self.Mi = ['%s.%d.tmp' %(self.f, i) for i in range(self.n_cpus)] # merged reads (temp)
        self.ci = ['q.%d.fst' %(i) for i in range(self.n_cpus)] # current reads
        self.Ci = ['q.%d.tmp' %(i) for i in range(self.n_cpus)] # current reads (temp)
        self.oi = ['otus.%d.fst' %(sid) for sid in self.sids] # otu representative sequences
        self.Oi = ['otus.%d.tmp' %(sid) for sid in self.sids] # otu representative sequences (temp)
        self.uc = ['otus.%d.uc' %(sid) for sid in self.sids] # uclust output files
        self.xi = ['otus.%d.counts' %(sid) for sid in self.sids] # otu tables (counts)
        
        # Get database for read mapping
        if self.denovo == True:
            self.db = self.oi
        elif self.ref_gg == True:
            self.db = ['%s/%d_otus.fasta' %(self.ggdb, sid) for sid in self.sids]
    
    def split_fastq(self):
        '''Split forward and reverse reads (for parallel processing)'''
        
        # Get list of commands
        cmds = []
        if self.f:
            cmd = 'python %s/split_fastq.py %s %s' %(self.library, self.f, self.n_cpus)
            cmds.append(cmd)
        if self.r:
            cmd = 'python %s/split_fastq.py %s %s' %(self.library, self.r, self.n_cpus)
            cmds.append(cmd)
        
        # Submit commands and validate output
        self.ssub.submit_and_wait(cmds, out=self.dry_run)
        self.ssub.validate_output(self.fi + self.ri, out=self.dry_run)
    
    def remove_primers(self):
        '''Remove diversity region + primer and discard reads with > 2 mismatches'''
        
        # Get list of commands
        cmds = []
        for i in range(self.n_cpus):
            if self.f:
                cmd = '%s/remove_primers.py %s %s --max_primer_diffs %d > %s' %(self.library, self.fi[i], self.p, self.p_mismatch, self.Fi[i])
                cmds.append(cmd)
            if self.r:
                cmd = '%s/remove_primers.py %s %s --max_primer_diffs %d > %s' %(self.library, self.ri[i], self.q, self.p_mismatch, self.Ri[i])
                cmds.append(cmd)
        
        # Submit commands and validate output
        self.ssub.submit_and_wait(cmds, out=self.dry_run)
        self.ssub.validate_output(self.Fi + self.Ri, out = self.dry_run)
        self.ssub.move_files(self.Fi + self.Ri, self.fi + self.ri, out = self.dry_run)
    
    def merge_reads(self):
        '''Merge forward and reverse reads using USEARCH'''
        
        # Intersect forward and reverse reads
        cmds = []
        for i in range(self.n_cpus):
            cmd = 'python %s/intersect_reads.py %s %s %s %s' %(self.library, self.fi[i], self.ri[i], self.Fi[i], self.Ri[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, out = self.dry_run)
        self.ssub.validate_output(self.Fi + self.Ri, out = self.dry_run)
        self.ssub.move_files(self.Fi + self.Ri, self.fi + self.ri, out = self.dry_run)
        
        # Merge reads
        cmds = []
        for i in range(self.n_cpus):
            cmd = '%s -fastq_mergepairs %s -reverse %s -fastq_truncqual %d -fastqout %s' %(self.usearch, self.fi[i], self.ri[i], self.truncqual, self.mi[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, out = self.dry_run)
        self.ssub.validate_output(self.mi, out = self.dry_run)
        self.ssub.remove_files(self.fi + self.ri, out = self.dry_run)
    
    def demultiplex_reads(self):
        '''Demultiplex samples using index and barcodes'''

        cmds = []
        for i in range(self.n_cpus):
            cmd = 'python %s/map_barcodes.py %s %s %s %d > %s' %(self.library, self.ci[i], self.b, self.x, self.b_mismatch, self.Ci[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        self.ssub.validate_output(self.Ci, self.dry_run)
        self.ssub.move_files(self.Ci, self.ci, self.dry_run)
    
    def quality_filter(self):
        '''Quality filter with truncqual and maximum expected error'''

        cmds = []
        for i in range(self.n_cpus):
            cmd = '%s -fastq_filter %s -fastq_truncqual %d -fastq_maxee %f -fastaout %s' %(self.usearch, self.ci[i], self.truncqual, self.maxee, self.Ci[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        self.ssub.validate_output(self.Ci, self.dry_run)
        self.ssub.move_files(self.Ci, self.ci, self.dry_run)
    
    def dereplicate(self):
        '''Concatenate files and dereplicate'''

        cmd = 'cat %s > q.fst' %(' '.join(self.ci))
        self.ssub.run_local([cmd], out = self.dry_run)
        self.ssub.validate_output(['q.fst'], self.dry_run)
        cmd = 'rm %s' %(' '.join(self.ci))
        self.ssub.run_local([cmd], out = self.dry_run)
        cmd = 'python %s/derep_fulllength.py q.fst q.derep.fst' %(self.library)
        self.ssub.submit_and_wait([cmd], self.dry_run)
        self.ssub.validate_output(['q.derep.fst'], self.dry_run)
    
    def denovo_clustering(rename=True):
        '''Denovo clustering with USEARCH'''

        cmds = []
        for i in range(len(self.sids)):
            sid = self.sids[i]
            cmd = '%s -cluster_otus q.derep.fst -otus %s -otuid .%d' %(self.usearch, self.oi[i], sid)
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        self.ssub.validate_output(self.oi, self.dry_run)
        
        # Rename OTUs
        if rename == True:
            cmds = []
            for i in range(len(self.sids)):
                sid = self.sids[i]
                cmd = 'python %s/usearch_python/fasta_number.py %s OTU%d_ > %s' %(self.library, self.oi[i], sid, self.Oi[i])
                cmds.append(cmd)
            self.ssub.submit_and_wait(cmds, self.dry_run)
            self.ssub.validate_output(self.Oi, self.dry_run)
            self.ssub.move_files(self.Oi, self.oi, self.dry_run)
    
    def remove_chimeras(self):
        '''Remove chimeras using gold database'''

        cmds = []
        for i in range(len(self.sids)):
            sid = self.sids[i]
            cmd = '%s -uchime_ref %s -db %s -nonchimeras %s -strand plus' %(self.usearch, self.oi[i], self.gold_db, self.Oi[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        self.ssub.validate_output(self.Oi, self.dry_run)
        self.ssub.move_files(self.Oi, self.oi, self.dry_run)
    
    def reference_mapping(self):
        '''Map reads to reference databases'''

        cmds = []
        for i in range(len(self.sids)):
            cmd = '%s -usearch_global q.derep.fst -db %s -uc %s -strand both -id .%d' %(self.usearch, self.db[i], self.uc[i], self.sids[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        self.ssub.validate_output(self.uc, self.dry_run)
    
    def make_otu_tables(self):
        '''Make OTU tables from UC file'''

        cmds = []
        for i in range(len(self.sids)):
            cmd = 'python %s/usearch_python/uc2otutab.py %s %s' %(self.library, self.uc[i], self.xi[i])
            cmds.append(cmd)
        self.ssub.submit_and_wait(cmds, self.dry_run)
        self.ssub.validate_output(self.xi, self.dry_run)
    

if __name__ == '__main__':
    # Initialize OTU caller
    oc = OTU_Caller()
    
    # Split fastq
    if not oc.dont_split:
        message('Splitting fastq')
        oc.split_fastq()
    
    # Remove primers
    if oc.primers == True:
        message('Removing primers')
        oc.remove_primers()
    
    # Merge reads
    if oc.merge == True:
        message('Merging reads')
        oc.merge_reads()
        
    # Set current reads
    if oc.merge == True:
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
        oc.dereplicate()
    
    # Denovo clustering
    if oc.denovo == True:
        message('Denovo clustering')
        oc.denovo_clustering(rename = True)
    
    # Map to reference database
    if oc.ref_gg == True:
        message('Mapping to reference')
        oc.reference_mapping()
    
    # Make OTU tables
    if oc.otu_table == True:
        message('Making OTU tables')
        oc.make_otu_tables()
