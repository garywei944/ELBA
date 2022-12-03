import sys
import fileinput as fi

def usage():
    sys.stderr.write("\nUsage: python {} <reads.fa>\n\n".format(sys.argv[0]))
    return -1

def main(argc, argv):

    if argc < 1: return usage()

    i, l = 1, 0
    name = ""
    for line in fi.input():
        if line.startswith(">"):
            if l > 0:
                sys.stdout.write("\t".join([str(i), str(l), name]) + "\n")
                i, l = i+1, 0
            name = line.lstrip(">").rstrip().split()[0]
        else:
            l += len(line.rstrip())
    if l > 0:
        sys.stdout.write("\t".join([str(i), str(l), name]) + "\n")

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
