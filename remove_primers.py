'''
Script that removes primers from reads in fastq files.

With a primer
    AAAAA

entries like
    @lolapolooza
    TAAAACATCATCATCAT
    +whatever
    abcdefghijklmnopq

become entries like
    @lolapolooza
    CATCATCATCAT
    +
    fghijklmnopq
'''

import argparse, itertools, sys
import util
#from util import *

class PrimerRemover():
    def __init__(self, fastq, primer, max_primer_diffs):
        '''
        Remove well-matched primers from sequences in a fastq file

        Parameters
        fastq : sequence or iterator of strings
            lines in the fastq file
        primer : string
            the primer sequence to be removed
        max_primer_diffs : int
            the maximum number of mismatches allowed before throwing out the read
        '''

        self.iterator = util.fastq_iterator(fastq)
        self.primer = primer
        self.primer_length = len(self.primer)
        self.max_primer_diffs = max_primer_diffs
        
        self.n_successes = 0
        self.n_failures = 0

    def entries(self):
        '''iterator over successfully trimmed input fastq entries'''

        for entry in self.iterator:
            [at_line, seq_line, quality_line] = self.iterator.next()

            # find the best primer position in the sequence
            primer_start_index, n_primer_diffs = util.mismatches(seq_line, self.primer, 15)

            # if we find a good match, trim the sequence and the
            # quality line and yield a single string
            if n_primer_diffs <= self.max_primer_diffs:
                primer_end_index = primer_start_index + self.primer_length
                trimmed_seq = seq_line[primer_end_index:]
                trimmed_quality = quality_line[primer_end_index:]
                self.n_successes += 1

                yield "\n".join([at_line, trimmed_seq, '+', trimmed_quality])
            else:
                self.n_failures += 1

    def print_entries(self):
        '''print the successfully trimmed entries'''

        timer = util.timer()
        for entry in self.entries():
            print entry
        self.elapsed_time = timer.next()


if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help='input fastq file')
    parser.add_argument('primer', help='primer sequence')
    parser.add_argument('-m', '--max_primer_diffs', default=0, type=int, help='maximum number of nucleotide mismatches in the primer (default: 0)')
    args = parser.parse_args()

    r = PrimerRemover(open(args.fastq), args.primer, args.max_primer_diffs)
    r.print_entries()