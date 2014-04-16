'''
split fastq into k files

usage:
  python split_fastq.py in.fsq 3

creates 3 files:
  in.fsq.1
  in.fsq.2
  in.fsq.3
  
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

for entry in iter_fsq(fst_fn):
    fh = fhs.next()
    fh.write('%s\n' %('\n'.join(entry)))