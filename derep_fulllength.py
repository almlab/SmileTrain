import sys

f = sys.argv[1]
o = sys.argv[2]
u = sys.argv[3]

out = open(o, 'w')
uc = open(u, 'w')

fst = {}

for line in open(f):
    line = line.rstrip()
    if line.startswith('>'):
        sid = line[1:]
    else:
        seq = line
        if seq not in fst:
            fst[seq] = [sid]
        else:
            fst[seq].append(sid)

i = 0
for seq in sorted(fst.keys(), key=lambda x: len(fst[x]), reverse=True):

    if len(fst[seq]) < 2:
        continue

    i += 1
    otu = 'seq%d' %(i)
    out.write( '>%s;size=%d\n%s\n' %(otu, len(fst[seq]), seq))
    uc.write('%s\t%s\n' %(otu, '\t'.join(fst[seq])))

out.close()
uc.close()
