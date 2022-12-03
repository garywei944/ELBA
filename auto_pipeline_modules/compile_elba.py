import sys
import getopt
import subprocess as sp
from pathlib import Path
import shutil
import os

def logmsg(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()

def abolsolutize(p):
    return p.absolute().resolve()

def usage():
    sys.stderr.write("Usage: python {} [options] /path/to/progdir\n\n".format(sys.argv[0]))
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -L INT,INT,...   prune k-mers below this bound [20]\n")
    sys.stderr.write("    -U INT,INT,...   prune k-mers above this bound [30]\n")
    sys.stderr.write("    -d FLOAT chernoff bound [0.1]\n")
    sys.stderr.write("    -p FILE  relative path to ELBA main directory [../]\n")
    return -1

def main(argc, argv):

    if argc < 2: return usage()

    elba_path = Path("..").absolute().resolve()
    lower_kmer_bounds = [20]
    upper_kmer_bounds = [30]
    delta_chernoff = 0.1

    try: opts, args = getopt.gnu_getopt(argv[1:], "L:U:d:p:h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        return usage()

    for o, a in opts:
        if o == "-h": return usage()
        elif o == "-L": lower_kmer_bounds = list(map(lambda x: int(x), a.split(",")))
        elif o == "-U": upper_kmer_bounds = list(map(lambda x: int(x), a.split(",")))
        elif o == "-d": delta_chernoff = float(a)

    if len(args) != 1:
        return usage()

    for l in L:
        for u in U:
            if l >= u:
                sys.stderr.write("error: can't have a lower k-mer bound greater than an upper k-mer bound\n")
                return usage()

    #outdir_path = absolutize(Path(args[0]))
    outdir_path = Path(args[0])
    outdir_path.mkdir(parents=True, exist_ok=True)

    logmsg("Created/touched output directory " + str(outdir_path))

    for p in outdir_path.iterdir():
        if p.is_dir():
            shutil.rmtree(str(p))
            logmsg("Removed directory '{}'".format(str(p)))
        elif p.is_file():
            p.unlink()
            logmsg("Removed file '{}'".format(str(p)))

    build_path = outdir_path.joinpath(
