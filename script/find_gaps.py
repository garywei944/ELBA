import sys
import numpy as np
from intervaltree import Interval, IntervalTree

LT_QUERY=0
EQ_QUERY=1
GT_QUERY=2

symbols = ["LT", "EQ", "GT"]

def get_contig_names(t, beg, end):
    intervals = sorted(t[beg:end])
    n = len(intervals)
    if n == 0: return "gap"
    else: return ",".join(it.data for it in intervals)

def main(depth_fname, query_depth, contigs_paf_fname):

    t = IntervalTree()
    for line in open(contigs_paf_fname, "r"):
        tokens = line.rstrip().split()
        contig_name = tokens[0]
        beg = int(tokens[7])
        end = int(tokens[8])
        t[beg:end] = contig_name

    length = 0
    pileup_coords = []

    for line in open(depth_fname, "r"):
        tokens = line.rstrip().split()
        pos = int(tokens[0])-1
        cov = int(tokens[1])
        length = max(length, pos+1)
        pileup_coords.append((pos, cov))

    pileup_vector = np.zeros(length, dtype=np.uint8) 

    for pos, cov in pileup_coords:
        pileup_vector[pos] = cov

    if pileup_vector[0] < query_depth: state = LT_QUERY
    elif pileup_vector[0] > query_depth: state = GT_QUERY
    else: state = EQ_QUERY

    sys.stdout.write("#relation\tstartpos\tlength\tcontig(s)\n")
    sys.stdout.flush()

    i = 0
    for j in range(1, length):
        if pileup_vector[j] < query_depth:
            if state != LT_QUERY:
                sys.stdout.write("{}\t{}\t{}\t{}\n".format(symbols[state], i+1, j-i, get_contig_names(t, i, j)))
                state = LT_QUERY
                i = j
        elif pileup_vector[j] > query_depth:
            if state != GT_QUERY:
                sys.stdout.write("{}\t{}\t{}\t{}\n".format(symbols[state], i+1, j-i, get_contig_names(t, i, j)))
                state = GT_QUERY
                i = j
        else:
            if state != EQ_QUERY:
                sys.stdout.write("{}\t{}\t{}\t{}\n".format(symbols[state], i+1, j-i, get_contig_names(t, i, j)))
                state = EQ_QUERY
                i = j

    sys.stdout.write("{}\t{}\t{}\n".format(symbols[state], i+1, j+1-i, get_contig_names(t, i, j+1)))
    sys.stdout.flush()

if __name__ == "__main__":

    if len(sys.argv) != 4:
        sys.stderr.write("usage: {} <depth.txt> <query depth> <contigs.paf>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    depth_fname = sys.argv[1]
    query_depth = int(sys.argv[2])
    contigs_paf_fname = sys.argv[3]

    main(depth_fname, query_depth, contigs_paf_fname)
