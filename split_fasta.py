'''
split fasta into k files

usage:
  python split_fasta.py in.fst 3

creates 3 files:
  in.fst.1
  in.fst.2
  in.fst.3
  
'''

import itertools, os.path, sys
from util import *

fst_fn = sys.argv[1]
k = int(sys.argv[2])

fns = ['%s.%d' %(fst_fn, i) for i in range(k)]

for fn in fns:
    if os.path.exists(fn):
        exit('file %s exists' %(fn))

fhs = cycle([open(fn, 'w') for fn in fns])

for [sid, seq] in iter_fst(fst_fn):
    fh = fhs.next()
    fh.write('>%s\n%s\n' %(sid, seq))