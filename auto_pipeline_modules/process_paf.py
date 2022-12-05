#!/usr/bin/env python

import sys
from pathlib import Path

def usage():
    sys.stderr.write("\nUsage: process_paf.py <in.paf> <idxtable.tsv>\n\n")
    return -1

def logmsg(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()

def fast_forward(d, f, name, length):
    n = ""
    while n != name:
        line = next(f)
        if not line:
            logmsg("error: idxtable file doesn't have entry for '{}'".format(name))
            sys.exit(-1)
        tokens = line.rstrip().split()
        old_name, l, new_name = tokens[0], int(tokens[1]), tokens[2]
        d[old_name] = (new_name, l)
        if old_name == name and l != length:
            logmsg("error: idxtable file has entry inconsistent with paf file for '{}'".format(name))
            sys.exit(-1)

def main(argc, argv):

    if argc < 3: return usage()

    paf_path = Path(argv[1]).resolve()
    idxtable_path = Path(argv[2]).resolve()

    if not paf_path.is_file():
        logmsg("error: {} is not a file".format(str(paf_path)))
        return usage()

    if not idxtable_path.is_file():
        logmsg("error: {} is not a file".format(str(idxtable_path)))
        return usage()

    pfh = open(str(paf_path), "r")
    ifh = open(str(idxtable_path), "r")

    d = {}

    for pafline in pfh:

        tokens = pafline.rstrip().split()

        query_name, target_name = tokens[0], tokens[5]
        query_length, target_length = int(tokens[1]), int(tokens[2])

        if not query_name in d:
            fast_forward(d, ifh, query_name, query_length)

        if not target_name in d:
            fast_forward(d, ifh, target_name, target_length)

        tokens[0] = d[query_name][0]
        tokens[5] = d[target_name][0]

        sys.stdout.write("\t".join(tokens) + "\n")

    sys.stdout.flush()

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
