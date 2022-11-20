import sys

def parse_paf_line(line):
    tokens = line.rstrip().split()
    name_a, name_b = tokens[0], tokens[5]
    len_a, beg_a, end_a = (int(v) for v in tokens[1:4])
    len_b, beg_b, end_b = (int(v) for v in tokens[6:9])
    return name_a, name_b, len_a, len_b, beg_a, beg_b, end_a, end_b

def parse_paf(paf_filename):

    read_lengths = {}
    read_intervals = {}

    def update(line, name, length, beg, end):
        if not name in read_lengths:
            read_lengths[name] = length
            read_intervals[name] = []
        read_intervals[name].append((beg, end))
        sys.stdout.flush()

    for line in open(paf_filename, "r"):
        name_a, name_b, len_a, len_b, beg_a, beg_b, end_a, end_b = parse_paf_line(line)
        update(line.rstrip(), name_a, len_a, beg_a, end_a)
        update(line.rstrip(), name_b, len_b, beg_b, end_b)

    return read_lengths, read_intervals

def find_chimera(paf_filename, coverage_min):

    read_lengths, read_intervals = parse_paf(paf_filename)

    for name, length in read_lengths.items():

        coverage = [0 for i in range(length)]

        for beg, end in read_intervals[name]:
            for i in range(beg, end):
                coverage[i] += 1

        in_gap = True

        middle_gaps = []
        extremity_gaps = []
        beg, end = 0, 0

        for i in range(len(coverage)):
            if coverage[i] <= coverage_min and not in_gap:
                beg, end = 0, 0
                beg = i
                in_gap = True

            if coverage[i] > coverage_min and in_gap:
                end = i
                in_gap = False
                if beg != end:
                    if beg == 0 or end == length:
                        extremity_gaps.append((beg, end))
                    else:
                        middle_gaps.append((beg, end))

        if in_gap:
            end = i + 1
            if beg != end:
                if beg == 0 or end == length:
                    extremity_gaps.append((beg, end))
                else:
                    middle_gaps.append((beg, end))

        if len(middle_gaps) > 0:
            sys.stdout.write("Chimeric:{},{};".format(name, length))
            for gap in middle_gaps:
                sys.stdout.write("{},{},{};".format(abs(gap[0]-gap[1]), gap[0], gap[1]))
            sys.stdout.write("\n")
            sys.stdout.flush()
        elif len(extremity_gaps) > 0:
            for gap in extremity_gaps:
                if abs(gap[0]-gap[1]) > 0.8 * length:
                    sys.stdout.write("Not_covered:{},{};".format(name, length))
                    sys.stdout.write("{},{},{};".format(abs(gap[0]-gap[1]), gap[0], gap[1]))
                    sys.stdout.write("\n")
                    sys.stdout.flush()
                    break;

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("usage: python yacrd.py [mapping.paf] [coverage_min]\n")
        sys.stderr.flush()
        sys.exit(1)

    paf_filename = sys.argv[1]
    coverage_min = int(sys.argv[2])

    find_chimera(paf_filename, coverage_min)


