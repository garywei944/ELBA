import sys
import getopt
import subprocess as sp
from pathlib import Path

def process_paf(in_paf_fname, out_paf_fname, idx_table_fname):

    in_paf = Path(in_paf_fname)
    out_paf = Path(out_paf_fname)

    idx_dict = dict(map(lambda line: line.rstrip().split(), open(idx_table_fname)))

    def nameize_pafline(pafline):
        l = pafline.rstrip().split()
        l[0], l[5] = idx_dict[l[0]], idx_dict[l[5]]
        return "\t".join(l) + "\n"

    with open(out_paf_fname, "w") as f:
        f.writelines(map(nameize_pafline, open(in_paf_fname)))

def write_idx_table(reads_fname, idx_table_fname):

    with open(idx_table_fname, "w") as f:
        p = sp.Popen(["awk", "/^>/{print ++i, substr($1, 2); next}", reads_fname], stdout=f)
        p.wait()

    return int(sp.Popen(["tail", "-n1", idx_table_fname], stdout=sp.PIPE).communicate()[0].decode().lstrip().rstrip().split()[0])

def usage():
    sys.stderr.write("Usage: python {} [options] <reads.fa> <elba>\n\n".format(sys.argv[0]))
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -t INT   ktip threshold [0]\n")
    sys.stderr.write("    -x INT   xdrop value [15]\n")
    sys.stderr.write("    -f STR   path prefixes [elba.out]\n")
    sys.stderr.write("    -n INT   number of processes [1]\n")
    sys.stderr.write("    -r FILE  reference fasta file\n")
    sys.stderr.write("    -p       prune bridges\n")
    sys.stderr.write("    -M       run on personal computer\n")
    return -1

def main(argc, argv):

    if argc < 3: return usage()

    path_prefix = "elba.out"
    prune_bridges = False
    ktip_threshold = 0
    num_procs = 1
    xdrop = 15
    reference_fname = ""
    on_mac = False

    try: opts, args = getopt.gnu_getopt(argv[1:], "pMt:f:n:x:h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        return usage()

    for o, a in opts:
        if o == "-h": return usage()
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
    contig_path = Path(path_prefix + ".ctg.fa").absolute().resolve()
    idx_table_path = Path(path_prefix + ".idxtable.txt").absolute().resolve()

    num_reads = write_idx_table(str(reads_path), str(idx_table_path))

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

    sys.stdout.write(" ".join(elba_cmd) + "\n")
    sys.stdout.flush()

    p = sp.Popen(elba_cmd)
    p.wait()

    for rankfile in Path.cwd().glob("elba_rank_*_log.txt"):
        rankfile.unlink()

    process_paf("elba.overlap.paf", path_prefix + ".overlap.paf", str(idx_table_path))
    process_paf("elba.string.paf", path_prefix + ".string.paf", str(idx_table_path))

    Path("elba.overlap.paf").unlink()
    Path("elba.string.paf").unlink()

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
