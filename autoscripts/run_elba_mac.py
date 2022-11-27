import sys
import getopt
import subprocess as sp
from pathlib import Path

def get_num_reads(reads_fname):
    grep_proc = sp.Popen(["grep", ">", reads_fname], stdout=sp.PIPE)
    wc_proc = sp.Popen(["wc", "-l"], stdin=grep_proc.stdout, stdout=sp.PIPE)
    return int(wc_proc.communicate()[0].decode().lstrip().rstrip())


def usage():
    sys.stderr.write("Usage: python {} [options] <reads.fa> <elba>\n\n".format(sys.argv[0]))
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -c FILE  contig fasta filename [ctg.fa]\n")
    sys.stderr.write("    -p       prune bridges\n")
    sys.stderr.write("    -t INT   ktip threshold [0]\n")
    sys.stderr.write("    -x INT   xdrop value [15]\n")
    sys.stderr.write("    -f STR   path prefixes [elba]\n")
    sys.stderr.write("    -n INT   number of processes [1]\n")
    return -1

def main(argc, argv):

    if argc < 3: return usage()

    contig_fname = "ctg.fa"
    path_prefix = "elba"
    prune_bridges = False
    ktip_threshold = 0
    num_procs = 1
    xdrop = 15

    try: opts, args = getopt.gnu_getopt(argv[1:], "c:pt:f:n:x:h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        return usage()

    for o, a in opts:
        if o == "-h": return usage()
        elif o == "-c": contig_fname = a
        elif o == "-p": prune_bridges = True
        elif o == "-t": ktip_threshold = int(a)
        elif o == "-f": path_prefix = a
        elif o == "-n": num_procs = int(a)
        elif o == "-x": xdrop = int(a)

    if len(args) != 2:
        return usage()

    reads_path = Path(args[0]).absolute().resolve()
    elba_exepath = Path(args[1]).absolute().resolve()
    contig_path = Path(contig_fname).absolute().resolve()

    num_reads = get_num_reads(str(reads_path))

    elba_cmd = ["mpirun", "-np", str(num_procs), str(elba_exepath), "-i", str(reads_path), "-o", str(contig_path), "--idxmap", "idmap", "-c", str(num_reads), "-k", "31", "--xa", str(xdrop), "-s", "1", "-O", "100000", "--afreq", "100000"]

    if prune_bridges:
        elba_cmd.append("--pb")

    if ktip_threshold != 0:
        elba_cmd += ["--tip", str(ktip_threshold)]

    output_fname = path_prefix + ".out"
    output_path = Path(output_fname).absolute().resolve()

    with open(str(output_path), "w") as f:
        p = sp.Popen(elba_cmd, stdout=f, stderr=sp.STDOUT)
        p.wait()


if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
