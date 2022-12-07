#!/usr/bin/env python

import sys
from pathlib import Path
import subprocess as sp

def usage():
    sys.stderr.write("\nUsage: process_paf.py <in.paf> <reads.fa>\n\n")
    return -1

def logmsg(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()

def fast_forward_closure(fidx):
    idx = 1
    def fast_forward(name_map, idx_name, length):
        nonlocal idx
        while not idx_name in name_map:
            try: line = next(fidx)
            except StopIteration:
                logmsg("error: found name in paf file that isn't in fasta file: '{}'".format(idx_name))
                sys.exit(usage())
            tokens = line.rstrip().split()
            read_name, read_length = tokens[0], int(tokens[1])
            name_map[str(idx)] = (read_name, read_length)
            idx += 1
        assert length == name_map[idx_name][1]
    return fast_forward

def main(argc, argv):

    if argc != 3: return usage()

    paf_path = Path(argv[1]).resolve()
    fasta_path = Path(argv[2]).resolve()

    for p in [paf_path, fasta_path]:
        if not p.is_file():
            logmsg("error: file '{}' not found".format(str(p)))
            return usage()

    proc = sp.Popen(["samtools", "faidx", str(fasta_path)])
    proc.wait()

    fasta_index_path = fasta_path.parent.resolve().joinpath(fasta_path.name + ".fai")

    if not fasta_index_path.is_file():
        logmsg("error: file '{}' should exist but doesn't".format(str(fasta_index_path)))
        return -1

    name_map = {}
    fidx = open(str(fasta_index_path), "r")


    ff = fast_forward_closure(fidx)

    for pafline in open(str(paf_path), "r"):

        tokens = pafline.rstrip().split()

        query_name, target_name = tokens[0], tokens[5]
        query_length, target_length = int(tokens[1]), int(tokens[6])

        ff(name_map, query_name, query_length)
        ff(name_map, target_name, target_length)

        tokens[0] = name_map[query_name][0]
        tokens[5] = name_map[target_name][0]

        sys.stdout.write("\t".join(tokens) + "\n")

    fidx.close()

    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))

