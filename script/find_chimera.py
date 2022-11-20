import sys
import numpy as np

def main(paf_fname, cov_min):

    pileup_vectors = {}

    for line in open(paf_fname, "r"):

        tokens = line.rstrip().split()
        query_name, target_name = tokens[0], tokens[5]
        query_length, target_length = int(tokens[1]), int(tokens[6])
        query_beg, target_beg = int(tokens[2]), int(tokens[7])
        query_end, target_end = int(tokens[3]), int(tokens[8])

        if not query_name in pileup_vectors:
            pileup_vectors[query_name] = np.zeros(query_length, dtype=np.uint8)

        if not target_name in pileup_vectors:
            pileup_vectors[target_name] = np.zeros(target_length, dtype=np.uint8)

        assert len(pileup_vectors[query_name]) == query_length
        assert len(pileup_vectors[target_name]) == target_length

        pileup_vectors[query_name][query_beg:query_end] += 1
        pileup_vectors[target_name][target_beg:target_end] += 1

    for read_name, pileup_vector in pileup_vectors.items():

        in_gap = True
        read_length = len(pileup_vector)
        gaps = []
        i = 0

        for j in range(read_length):

            if pileup_vector[j] <= cov_min and not in_gap:
                i = j
                in_gap = True

            if pileup_vector[j] > cov_min and in_gap:
                in_gap = False
                if i != j and i != 0:
                    gaps.append((i, j))

        if len(gaps) > 0:
            sys.stdout.write("Chimeric\t{}\t{}".format(read_name, read_length))
            for i,j in gaps:
                sys.stdout.write("\t{},{},{}".format(abs(i-j), i, j))
            sys.stdout.write("\n")
            sys.stdout.flush()

if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.stderr.write("usage: {} <map.paf> <coverage_min>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    paf_fname = sys.argv[1]
    cov_min = int(sys.argv[2])

    main(paf_fname, cov_min)
