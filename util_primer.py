import re, string, sys, time, itertools, os, subprocess
from Bio import SeqRecord, SeqIO, Seq
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
    def __init__(self, fastq, primer, max_primer_diffs, reverse_primer=None, skip=1, out=sys.stdout):
        '''
        Remove well-matched primers from sequences in a fastq file. If given a reverse
        primer, search from the beginning for the forward primer and from the end for
        the revese primer.

        Parameters
        fastq : filename or filehandle
            input
        primer : string
            the primer sequence to be removed
        max_primer_diffs : int
            the maximum number of mismatches allowed before throwing out the read
        reverse_primer : string
            a primer to be trimmed from the end of the sequence
        skip : integer >= 1
            take only every n-th entry (so skip=1 means every entry)
        '''

        self.records = SeqIO.parse(fastq, 'fastq')
        self.primer = primer
        self.primer_length = len(self.primer)
        self.max_primer_diffs = max_primer_diffs

        if reverse_primer is not None:
            #self.reverse_primer = Seq.Seq(reverse_primer).reverse_complement()
            self.reverse_primer = Seq.Seq(reverse_primer)
            self.reverse_primer_length = len(self.reverse_primer)
        else:
            self.reverse_primer = None

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
                primer_start_index, n_primer_diffs = mismatches(record.seq, self.primer, 15)

                quality = record.letter_annotations['phred_quality']
    
                # if we find a good match, trim the sequence and the
                # quality line and yield a single string
                if n_primer_diffs <= self.max_primer_diffs:
                    primer_end_index = primer_start_index + self.primer_length
                    trim_seq = record.seq[primer_end_index:]
                    trim_quality = quality[primer_end_index:]
                    success = True
                else:
                    success = False

                if success and self.reverse_primer is not None:
                    rc_trim_seq = trim_seq.reverse_complement()
                    reverse_start_index, n_reverse_diffs = mismatches(rc_trim_seq, self.reverse_primer, 15)

                    if n_reverse_diffs <= self.max_primer_diffs:
                        reverse_end_index = reverse_start_index + self.reverse_primer_length
                        rc_trim_seq = rc_trim_seq[reverse_end_index:]
                        trim_seq = rc_trim_seq.reverse_complement()
                        trim_quality = trim_quality[:-reverse_end_index]
                        success = True
                    else:
                        success = False

                if success:
                    new_record = SeqRecord.SeqRecord(seq=str(trim_seq), id=record.id, letter_annotations={'phred_quality': trim_quality}, description="")
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