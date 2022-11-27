import sys
import numpy as np

LT_QUERY=0
EQ_QUERY=1
GT_QUERY=2

symbols = ["<", "=", ">"]

def main(depth_fname, query_depth):

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

    i = 0
    for j in range(1, length):
        if pileup_vector[j] < query_depth:
            if state != LT_QUERY:
                sys.stdout.write("{}\t{}\t{}\n".format(i+1, j-i, symbols[state]))
                state = LT_QUERY
                i = j
        elif pileup_vector[j] > query_depth:
            if state != GT_QUERY:
                sys.stdout.write("{}\t{}\t{}\n".format(i+1, j-i, symbols[state]))
                state = GT_QUERY
                i = j
        else:
            if state != EQ_QUERY:
                sys.stdout.write("{}\t{}\t{}\n".format(i+1, j-i, symbols[state]))
                state = EQ_QUERY
                i = j

    sys.stdout.write("{}\t{}\t{}\n".format(i+1, j+1-i, symbols[state]))


if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.stderr.write("usage: {} <depth.txt> <query depth>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    depth_fname = sys.argv[1]
    query_depth = int(sys.argv[2])

    main(depth_fname, query_depth)
