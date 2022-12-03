import sys
import getopt
from pathlib import Path
import os

def usage():
    sys.stderr.write("\nUsage: python {} [options] <reads.fa> /path/to/elba\n\n".format(sys.argv[0]))
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -p       prune bridges\n")
    sys.stderr.write("    -t INT   ktip threshold [0]\n")
    sys.stderr.write("    -f STR   output file prefixes [elba]\n")
    sys.stderr.write("    -n INT   number of processes [1]\n\n")
    return -1

def main(argc, argv):

    if argc < 3: return usage()

    file_prefix = "elba"

    prune_bridges = False
    ktip_threshold = 0
    num_procs = 1

    try: opts, args = getopt.gnu_getopt(argv[1:], "pt:f:n:h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        return usage()

    for o, a in opts:
        if o == "-h": return usage()
        elif o == "-p": prune_bridges = True
        elif o == "-t": ktip_threshold = int(a)
        elif o == "-f": file_prefix = a
        elif o == "-n": num_procs = int(a)

    if len(args) != 2:
        return usage()

    reads_fname = Path(args[0]).resolve()
    elba_pname = Path(args[1]).resolve()




if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
