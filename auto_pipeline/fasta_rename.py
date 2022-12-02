import sys
import getopt

line_width_default = 0

def usage():

    sys.stderr.write("Usage: python {} [options] <in.fa>\n".format(sys.argv[0]))
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -D FILE   dump read name map\n")
    sys.stderr.write("    -l INT    number of residues per line; 0 for entire sequence [{}]\n".format(line_width_default))
    sys.stderr.write("    -n        include original read name as comment\n")
    sys.stderr.write("    -0        use zero-based indexing\n")
    sys.stderr.write("    -h        print help message\n")
    sys.stderr.flush()

    return -1


def main(argc, argv):

    if argc < 2: return usage()

    readmap_fname = None
    line_width = line_width_default
    include_comment = False
    one_based_indexing = True

    try: opts, args = getopt.gnu_getopt(argv[1:], "D:l:n0h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        sys.stderr.flush()
        return usage()

    for o, a in opts:
        if o == "-h": return usage()
        elif o == "-D": readmap_fname = a
        elif o == "-l": line_width = max(0, int(a))
        elif o == "-n": include_comment = True
        elif o == "-0": one_based_indexing = False

    if len(args) == 0:
        sys.stderr.write("error: missing <in.fa> [fatal]\n")
        sys.stderr.flush()
        return usage()

    fasta_fname = args[0]

    readmap = None if readmap_fname is None else open(readmap_fname, "w")
    idx = int(one_based_indexing)
    seqlines = []
    name_orig = ""

    def write_fasta_record():
        nonlocal seqlines, name_orig, readmap, idx, line_width, include_comment
        seq = "".join(seqlines)
        name = name_orig.split(None, 1)[0]
        if not readmap is None: readmap.write(name + "\n")
        if include_comment: sys.stdout.write(">{}\t{}\n".format(idx, name))
        else: sys.stdout.write(">{}\n".format(idx))
        idx += 1
        if line_width > 0:
            n = len(seq)
            for i in range(0, n, line_width): sys.stdout.write("{}\n".format(seq[i:min(i+line_width,n)]))
        else: sys.stdout.write("{}\n".format(seq))

    for line in open(fasta_fname, "r"):
        sys.stderr.write(str(idx) + "\n")
        if line.startswith(">"):
            if len(seqlines) > 0:
                write_fasta_record()
                seqlines = []
            name_orig = line.lstrip(">").rstrip()
        else: seqlines.append(line.rstrip())
    if len(seqlines) > 0: write_fasta_record()

    if not readmap is None: readmap.close()

    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
