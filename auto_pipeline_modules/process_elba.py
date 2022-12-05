#!/usr/bin/env python

import sys
import getopt
from pathlib import Path
import subprocess as sp
import os

def logmsg(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()

def usage():
    sys.stderr.write("\nUsage: process_elba.py [options] /path/to/elba_output\n\n")
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -r FILE  reference fasta for quast assembly metrics [null]\n")
    sys.stderr.write("    -q FILE  compile important quast metrics into this file [null]\n\n")
    return -1

def main(argc, argv):

    if argc < 2: return usage()

    reference_fname = ""
    quast_report_fname = ""

    try: opts, args = getopt.gnu_getopt(argv[1:], "r:q:h")
    except getopt.GetoptError as err:
        logmsg("error: {}".format(err))
        return usage()

    for o, a in opts:
        if o == "-h": return usage()
        if o == "-r": reference_fname = a
        if o == "-q": quast_report_fname = a

    outdir_path = Path(args[0]).resolve()

    if not outdir_path.is_dir():
        logmsg("error: there is no directory at '{}'".format(str(outdir_path)))
        return usage()

    for logfile in outdir_path.glob("elba_rank_*.txt"):
        logfile.unlink()


    if reference_fname != "":
        reference_path = Path(reference_fname).resolve()
        if not reference_path.is_file():
            logmsg("error: file not found at '{}'".format(str(reference_path)))
            return usage()

        quast_outdir = outdir_path.joinpath("quast_outdir").resolve()

        matches = list(outdir_path.glob("*.contigs.fa"))

        if len(matches) != 1:
            logmsg("error: couldn't find a fasta of contigs in '{}'".format(str(outdir_path)))
            return -1

        contigs_path = matches[0].resolve()

        if not contigs_path.is_file():
            logmsg("error: file not found at '{}'".format(str(contigs_path)))
            return -1

        quast_cmd = ["quast.py", "-o", str(quast_outdir), "-r", str(reference_path), str(contigs_path)]
        logmsg("executing: '{}'\n".format(" ".join(quast_cmd)))
        proc = sp.Popen(quast_cmd)
        proc.wait()


if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
