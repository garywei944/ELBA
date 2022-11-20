import sys

def read_fasta_name_map(fasta_fname):
    read_name_map = []
    for line in open(fasta_fname, "r"):
        if line.startswith(">"):
            name = line.lstrip(">").rstrip().split(None, 1)[0]
            read_name_map.append(name)
    return read_name_map

def main(fasta_fname, paf_fname):

    read_name_map = read_fasta_name_map(fasta_fname)

    for line in open(paf_fname, "r"):
        items = line.rstrip().split()
        items[0] = read_name_map[int(items[0])-1]
        items[5] = read_name_map[int(items[5])-1]
        sys.stdout.write("\t".join(items) + "\n")

    sys.stdout.flush()

    return 0

if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.stderr.write("usage: {} <reads.fa> <elba.paf>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    fasta_fname = sys.argv[1]
    paf_fname = sys.argv[2]

    retval = main(fasta_fname, paf_fname)

    sys.exit(retval)
