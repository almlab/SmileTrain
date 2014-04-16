import sys

f = sys.argv[1]
r = sys.argv[2]

of = sys.argv[3]
rf = sys.argv[4]

q = {}
z = {}

i=0
for line in open(f):
    line = line.rstrip()
    i += 1
    if i%4 == 1:
        sid = line[1:].split()[0]
        q[sid] = ''

i=0
for line in open(r):
    line = line.rstrip()
    i += 1
    if i%4 == 1:
        sid = line[1:].split()[0]
        if sid in q:
            z[sid] = 1

fout = open(of, 'w')
i=0
for line in open(f):
    line = line.rstrip()
    i += 1
    if i%4 == 1:
        sid = line[1:]
        kid = sid.split()[0]
    elif i%4 == 2:
        seq = line
    elif i%4 == 0:
        qua = line
        if kid in z:
            fout.write( '@%s\n%s\n+\n%s\n' %(sid, seq, qua) )

rout = open(rf, 'w')
i=0
for line in open(r):
    line = line.rstrip()
    i += 1
    if i%4 == 1:
        sid = line[1:]
        kid = sid.split()[0]
    elif i%4 == 2:
        seq = line
    elif i%4 == 0:
        qua = line
        if kid in z:
            rout.write( '@%s\n%s\n+\n%s\n' %(sid, seq, qua) )
