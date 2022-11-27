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
    sys.stderr.write("    -t INT   ktip threshold [0]\n")
    sys.stderr.write("    -x INT   xdrop value [15]\n")
    sys.stderr.write("    -f STR   path prefixes [elba]\n")
    sys.stderr.write("    -n INT   number of processes [1]\n")
    sys.stderr.write("    -p       prune bridges\n")
    sys.stderr.write("    -M       run on personal computer\n")
    return -1

def main(argc, argv):

    if argc < 3: return usage()

    contig_fname = "ctg.fa"
    path_prefix = "elba"
    prune_bridges = False
    ktip_threshold = 0
    num_procs = 1
    xdrop = 15
    on_mac = False

    try: opts, args = getopt.gnu_getopt(argv[1:], "c:pMt:f:n:x:h")
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
        elif o == "-M": on_mac = True

    if len(args) != 2:
        return usage()

    reads_path = Path(args[0]).absolute().resolve()
    elba_exepath = Path(args[1]).absolute().resolve()
    contig_path = Path(contig_fname).absolute().resolve()

    num_reads = get_num_reads(str(reads_path))

    elba_cmd = []

    if on_mac: elba_cmd += ["mpirun", "-np"]
    else: elba_cmd += ["srun", "-n"]

    elba_cmd.append(str(num_procs))

    if not on_mac:
        elba_cmd += ["-c", str(256//num_procs), "--cpu_bind=cores"]

    elba_cmd += [str(elba_exepath), "-i", str(reads_path), "-o", str(contig_path), "--idxmap", "idmap", "-c", str(num_reads), "-k", "31", "--xa", str(xdrop), "-s", "1", "-O", "100000", "--afreq", "100000"]

    if prune_bridges:
        elba_cmd.append("--pb")

    if ktip_threshold != 0:
        elba_cmd += ["--tip", str(ktip_threshold)]

    #  output_fname = path_prefix + ".out"
    #  output_path = Path(output_fname).absolute().resolve()

    sys.stdout.write(" ".join(elba_cmd) + "\n")
    sys.stdout.flush()

    p = sp.Popen(elba_cmd)
    p.wait()

    for rankfile in Path.cwd().glob("elba_rank_*_log.txt"):
        rankfile.unlink()

    Path("elba.overlap.paf").rename(Path(path_prefix + ".overlap.paf"))
    Path("elba.string.paf").rename(Path(path_prefix + ".string.paf"))

    #with open(str(output_path), "w") as f:
    #    p = sp.Popen(elba_cmd, stdout=f, stderr=sp.STDOUT)
    #    p.wait()

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
