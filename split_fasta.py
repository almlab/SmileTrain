'''
Split an input fasta file into k files.

usage:
  python split_fasta.py in.fst 3

creates 3 files:
  in.fst.1
  in.fst.2
  in.fst.3
  
'''

import itertools, os.path, sys, argparse
import util

def output_filenames(input_filename, k):
    '''destination filenames foo.fastq.0, etc.'''
    return ['%s.%d' % (input_filename, i) for i in range(k)]

def split_fasta_entries(fasta, filenames, by_hash=False):
    '''
    Send entries in the input to filenames, cycling over each filename.
    
    fasta : list or iterator of strings
        lines in the input fasta file
    filenames : list or iterator of strings
        output filenames
    by_hash : bool (default false)
        split based on hash value rather than just cycling
        
    returns : nothing
    '''
    
    # open all the filehandles
    fhs = [open(filename, 'w') for filename in filenames]
    
    if by_hash:
        fastas = util.fasta_entries(fasta, output_type='list')
        
        # make the function for determining which bin the sequences fall in
        h = lambda seq: hash(seq) % len(fhs)
        
        # pick the filehandle bashed on the sequence's hash, then write
        for sid, seq in fastas:
            fhs[h(seq)].write(">%s\n%s\n" %(sid, seq))    
    else:
        fh_cycler = util.cycle(fhs)
        
        # prepare an iterator over the fastq entries
        fastas = util.fasta_entries(fasta, output_type='string')
        
        for entry, fh in itertools.izip(fastas, fh_cycler):
            fh.write("%s\n" % entry)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split a fasta file foo.fasta into multiple fastq files foo.fasta.0, foo.fasta.1, etc.')
    parser.add_argument('fasta', help='input fasta')
    parser.add_argument('n_files', type=int, help='number of split files to output')
    parser.add_argument('-s', '--hash', action='store_true', help='split by hash')
    args = parser.parse_args()
    
    filenames = output_filenames(args.fasta, args.n_files)
    
    util.check_for_collisions(filenames)
    
    if len(filenames) == 1:
        # just copy the file
        shutil.copy(args.fasta, filenames[0])
    else:
        # split the file entry by entry
        with open(args.fasta) as f:
            split_fasta_entries(f, filenames, by_hash=args.hash)