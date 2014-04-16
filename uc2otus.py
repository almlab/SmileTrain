# Convert .uc file to OTU table

def parse_line(line):
    # Parse otu id and sample id from .uc line
    line = line.rstrip().split()
    otu = line[-1]
    sid = line[-3].split('_')[0]
    return [otu, sid]

def map_file(fn):
    # Convert .uc file to counts
    counts = {}
    sids = []
    otus = []
    for line in open(fn):
        if line.startswith('H'):
            [otu, sid] = parse_line(line)
            if otu and otu not in otus:
                otus.append(otu)
            if sid and sid not in sids:
                sids.append(sid)
            if otu and sid:
                if sid not in counts:
                    counts[sid] = {}
                if otu not in counts[sid]:
                    counts[sid][otu] = 0
                counts[sid][otu] += 1
    return counts, sids, otus

def print_counts(counts, sids, otus):
    # Print counts data to stdout
    sids = sorted(sids)
    otus = sorted(otus)
    print '\t' + '\t'.join(otus)
    for sid in sids:
        outline = [sid]
        for otu in otus:
            if otu in counts[sid]:
                outline.append('%d' %(counts[sid][otu]))
            else:
                outline.append('0')
        print '\t'.join(outline)

def run():
    import sys
    fn = sys.argv[1]
    counts, sids, otus = map_file(fn)
    print_counts(counts, sids, otus)

run()
