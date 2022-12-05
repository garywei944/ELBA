#!/usr/bin/env python

import sys
import getopt
from pathlib import Path
import subprocess as sp
import os

def count_reads(fname):
    grep_proc = sp.Popen(["grep", ">", fname], stdout=sp.PIPE)
    wc_proc = sp.Popen(["wc", "-l"], stdin=grep_proc.stdout, stdout=sp.PIPE)
    count = int(wc_proc.communicate()[0].decode())
    return count

def usage():
    sys.stderr.write("\nUsage: run_elba.py [options] <reads.fa> /path/to/elba\n\n")
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -p       prune bridges\n")
    sys.stderr.write("    -k INT   k-mer size [31]\n")
    sys.stderr.write("    -l INT   ktip threshold [0]\n")
    sys.stderr.write("    -x INT   xdrop value [15]\n\n")

    sys.stderr.write("    -o STR   output file prefixes [elba]\n")
    sys.stderr.write("    -J STR   job name [elba.saul]\n\n")

    sys.stderr.write("    -t INT   time limit in minues [30]\n")
    sys.stderr.write("    -n INT   number of processes [1]\n")
    sys.stderr.write("    -N INT   number of compute nodes [1] (overrides -n)\n\n")

    sys.stderr.write("    -u STR   email for user updates []\n\n")
    return -1

def perlmutter_task_resources(num_procs, num_nodes):

    proc_counts = [p**2 for p in range(1,12)] + [p**2 for p in range(1,512) if (p**2)%128==0]
    node_counts = [max(P//128, 1) for P in proc_counts]

    if num_nodes > 1:
        for i in range(len(node_counts)):
            if num_nodes <= node_counts[i]:
                break
        return proc_counts[i], node_counts[i], 2
    else:
        for i in range(len(proc_counts)):
            if num_procs <= proc_counts[i]:
                break
        return proc_counts[i], node_counts[i], max(256//proc_counts[i], 2)

def main(argc, argv):

    if argc < 3: return usage()

    file_prefix = "elba"
    jobname = "elba.saul"
    email = ""
    prune_bridges = False
    ktip_threshold = 0
    time_limit = 30
    num_procs = 1
    num_nodes = 1
    kmer_size = 31
    xdrop = 15

    try: opts, args = getopt.gnu_getopt(argv[1:], "pl:t:o:n:N:J:u:x:k:h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        return usage()

    for o, a in opts:
        if o == "-h": return usage()
        elif o == "-p": prune_bridges = True
        elif o == "-l": ktip_threshold = int(a)
        elif o == "-t": time_limit = int(a)
        elif o == "-o": file_prefix = a
        elif o == "-J": jobname = a
        elif o == "-u": email = a
        elif o == "-n": num_procs = int(a)
        elif o == "-N": num_nodes = int(a)
        elif o == "-k": kmer_size = int(a)
        elif o == "-x": xdrop = int(a)

    if len(args) != 2:
        return usage()

    reads_fname = Path(args[0]).resolve()
    elba_pname = Path(args[1]).resolve()

    num_reads = count_reads(str(reads_fname))

    num_procs, num_nodes, cpus_per_proc = perlmutter_task_resources(num_procs, num_nodes)

    qos = "debug" if time_limit <= 30 and num_nodes == 1 else "regular"

    sys.stdout.write("#!/bin/bash\n\n")
    sys.stdout.write("#SBATCH -N {}\n".format(num_nodes))
    sys.stdout.write("#SBATCH -C cpu\n")
    sys.stdout.write("#SBATCH -q {}\n".format(qos))
    sys.stdout.write("#SBATCH -J {}\n".format(jobname))
    sys.stdout.write("#SBATCH -t {}\n".format(time_limit))
    sys.stdout.write("#SBATCH --error={}.%j.err\n".format(jobname))
    sys.stdout.write("#SBATCH --output={}.%j.out\n".format(jobname))
    sys.stdout.write("#SBATCH --switches=1\n")

    if email != "":
        sys.stdout.write("#SBATCH --mail-user={}\n".format(email))
        sys.stdout.write("#SBATCH --mail-type=ALL\n")

    sys.stdout.write("\nexport OMP_NUM_THREADS=1\n")
    sys.stdout.write("export OMP_PLACES=threads\n")
    sys.stdout.write("export OMP_PROC_BIND=spread\n\n")

    cmd = ["srun",
           "-n", str(num_procs),
           "-N", str(num_nodes),
           "-c", str(cpus_per_proc),
           "--cpu_bind=cores",
           str(elba_pname),
           "-i", str(reads_fname),
           "-o", file_prefix,
           "-k", str(kmer_size),
           "--idxmap", "idmap",
           "-c", str(num_reads),
           "--af", "{}.af".format(jobname),
           "--xa", str(xdrop),
           "-s", str(1),
           "-O", str(100000),
           "--afreq", str(100000)]

    if prune_bridges: cmd += ["--pb"]

    if ktip_threshold > 0: cmd += ["--tip", str(ktip_threshold)]


    sys.stdout.write(" ".join(cmd) + "\n")
    sys.stdout.flush()

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
