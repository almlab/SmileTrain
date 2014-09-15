import re, string, sys, time, itertools, os, subprocess
from Bio import SeqRecord, SeqIO
from SmileTrain import util
import usearch_python.primer

def mismatches(seq, primer, w):
    '''
    Calculate mismatches between a sequence and primer with window size w.
    Returns the starting index and number of mismatches for the best match.
    
    Parameters
    seq : string (or stringable)
        target (longer) sequence
    primer : string (or stringable)
        query (shorter) sequence
    w : int
        window size (how many bases on the target to go in, at most)
        
    returns : tuple
        (primer start index, differences in primer and sequence)
    '''

    I = 0
    D = len(seq)
    for i in range(w):
        d = usearch_python.primer.MatchPrefix(seq[i:], primer)
        if d < D:
            I = i
            D = d
    return (I, D)

class PrimerRemover():
    def __init__(self, fastq, primer, max_primer_diffs, output_type='string', skip=1, out=sys.stdout):
        '''
        Remove well-matched primers from sequences in a fastq file

        Parameters
        fastq : filename or filehandle
            input
        primer : string
            the primer sequence to be removed
        max_primer_diffs : int
            the maximum number of mismatches allowed before throwing out the read
        output_type : string 'string' or 'list' (default 'string')
            output format for iterator. string is a single string; list is the 3-element list
        skip : integer >= 1
            take only every n-th entry (so skip=1 means every entry)
        '''

        self.records = SeqIO.parse(fastq, 'fastq')
        self.primer = primer
        self.primer_length = len(self.primer)
        self.max_primer_diffs = max_primer_diffs

        self.output_type = output_type
        self.skip = skip
        self.out = out
        
        self.n_successes = 0
        self.n_failures = 0
        self.position = 0

    def __iter__(self):
        return self

    def next(self):
        '''iterator over successfully trimmed input fastq entries'''

        while True:
            record = self.records.next()

            # take this entry if it's the n-th entry
            if self.position % self.skip == 0:
                # find the best primer position in the sequence
                primer_start_index, n_primer_diffs = mismatches(str(record.seq), self.primer, 15)
    
                # if we find a good match, trim the sequence and the
                # quality line and yield a single string
                if n_primer_diffs <= self.max_primer_diffs:
                    primer_end_index = primer_start_index + self.primer_length
                    new_record = SeqRecord.SeqRecord(seq=str(record.seq)[primer_end_index:], id=record.id, letter_annotations={'phred_quality': record.letter_annotations['phred_quality'][primer_end_index:]})
                    self.n_successes += 1
                    
                    return new_record
                else:
                    self.n_failures += 1
                
            self.position += 1
            
    def print_entries(self):
        '''print the successfully trimmed entries'''

        timer = util.timer()
        for record in self:
            SeqIO.write(record, self.out, 'fastq')

        self.elapsed_time = timer.next()
        
    def check_entries(self):
        '''check for successfully trimmed entries'''
        
        timer = util.timer()
        for entry in self: pass
        self.elapsed_time = timer.next()

    def diagnostic_message(self):
        '''return a string about how the trimming went'''
        success_rate = float(self.n_successes) / (self.n_successes + self.n_failures) * 100
        return "Primer removal:\n%s successes and %s failures (%.0f%% success rate)\n%f seconds elapsed\n" %(self.n_successes, self.n_failures, success_rate, self.elapsed_time)