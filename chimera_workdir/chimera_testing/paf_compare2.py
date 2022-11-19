import sys

def read_paf(paf_filename):
    overlaps = {}
    for line in open(paf_filename, "r"):
        items = line.rstrip().split()
        rowid = int(items[0])
        rowlen = int(items[1])
        rowbeg = int(items[2])
        rowend = int(items[3])
        strand = items[4]
        colid = int(items[5])
        collen = int(items[6])
        colbeg = int(items[7])
        colend = int(items[8])

        if not rowid in overlaps:
            overlaps[rowid] = {}

        overlaps[rowid][colid] = (line.rstrip(), rowid, rowlen, rowbeg, rowend, strand, colid, collen, colbeg, colend)

    return overlaps

mm_overlaps = read_paf("mm.paf")
elba_overlaps = read_paf("elba.paf")

for rowid in mm_overlaps:
    if rowid in elba_overlaps:
        for colid in mm_overlaps[rowid]:
            if colid in elba_overlaps[rowid]:
                o1 = elba_overlaps[rowid][colid]
                o2 = mm_overlaps[rowid][colid]
                sys.stdout.write("{}\t{}\t::\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}\n".format(rowid,
                                                                                     colid,
                                                                                     o1[2]-o2[2],
                                                                                     o1[3]-o2[3],
                                                                                     o1[4]-o2[4],
                                                                                     o1[5],
                                                                                     o2[5],
                                                                                     o1[7]-o2[7],
                                                                                     o1[8]-o2[8],
                                                                                     o1[9]-o2[9]))
                #  sys.stdout.write("ELBA\t{}\n".format(elba_overlaps[rowid][colid][0]))
                #  sys.stdout.write("MMAP\t{}\n.\n".format(mm_overlaps[rowid][colid][0]))


