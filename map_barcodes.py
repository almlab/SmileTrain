import primer
import sys
from string import maketrans

rctab = maketrans('ACGTacgt','TGCAtgca')

def reverse_complement(x):
    return x[::-1].translate(rctab)

def mismatches(seq, subseq, w):
    I = 0
    D = len(seq)
    for i in range(w):
        d = primer.MatchPrefix(seq[i:], subseq)
        if d < D:
            I = i
            D = d
    return [I, D]

# Read input arguments
fasta = sys.argv[1]
bcode = sys.argv[2]
index = sys.argv[3]
MAX_BARCODE_DIFFS = int(sys.argv[4])

# Map each barcode to a sample
b2s = {}
for line in open(bcode):
    s, b = line.rstrip().split()
    b2s[reverse_complement(b)] = s

# Get list of reads
reads = {}
k = 0
for line in open(fasta):
    line = line.rstrip()
    k += 1
    if k%4 == 1:
        reads[line[1:]] = {}

# Map each read to a sample
# (read -> index -> barcode -> sample)
r2s = {}
k = 0
for line in open(index):
    line = line.rstrip()
    k += 1
    if k%4 == 1:
        r = line[1:]
        if r not in reads:
            r = ''
    elif r != '' and k%4 == 2:
        s = line

        # map to best barcode
        B = ''
        D = ''
        for b in b2s.keys():
            q, d = mismatches(s, b, 1)
            if d < D:
                B = b
                D = d
        if D <= MAX_BARCODE_DIFFS:
            r2s[r] = b2s[B]

count = {}

success = 0
fail = 0

for line in open(fasta):
    line = line.rstrip()
    if line.startswith('>'):
        read = line[1:]
        if read not in r2s:
            read = ''
            fail += 1
            continue
        sample = r2s[read]
        if sample not in count:
            count[sample] = 0
        count[sample] += 1
    elif read != '':
        print '>barcode=%s;%d' %(sample, count[sample])
        print line
        success += 1
