#!/usr/bin/env python

import sys
import getopt
import subprocess as sp
from pathlib import Path
import shutil
import os
import platform

command_id = 1

def logmsg(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()

def usage():
    sys.stderr.write("\nUsage: compile_elba.py [options] /path/to/progdir\n\n")
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -L INT,INT,...   prune k-mers below this bound [20]\n")
    sys.stderr.write("    -U INT,INT,...   prune k-mers above this bound [30]\n")
    sys.stderr.write("    -d FLOAT         chernoff bound [0.1]\n")
    sys.stderr.write("    -I PATH          compilation directory [compdir]\n")
    sys.stderr.write("    -h               give more detail in example form\n\n")
    return -1

def usage_extra():
    usage()
    sys.stderr.write("Example:\n")
    sys.stderr.write("    python {} -L10,15,20 -U40,35,30 -d0.2 -I compdir progdir\n\n".format(sys.argv[0]))
    sys.stderr.write("    This sets up compilation files in $compdir and outputs 3 executables to $progdir\n"
                     "    named elba.l10.u40, elba.l15.u35, and elba.l20,u30. These are ELBA executables with\n"
                     "    lower and upper k-mer bounds coming after the l and the u respectively. Delta chernoff\n"
                     "    value is set to 0.2 on all of them\n\n")
    return -1

def run_command(command_list):
    global command_id
    step_id = 1
    logmsg("\n({}.*) executing: '{}'\n".format(command_id, " ".join(command_list)))
    proc = sp.Popen(command_list, stdout=sp.PIPE)
    while True:
        output = proc.stdout.readline().decode()
        if output == "" and proc.poll() is not None:
            break
        if output:
            logmsg("\t({}.{}). {}".format(command_id, step_id, output.rstrip()))
            step_id += 1
    command_id += 1

def main(argc, argv):

    osplat = platform.system()

    if not osplat in {"Linux", "Darwin"}:
        logmsg("error: '{}' is not a valid platform to run this script on. Must be Linux or Darwin".format(osplat))
        return -1

    envnames = ["ELBA_HOME", "COMBBLAS_HOME"]
    envpaths = []

    for i, envname in enumerate(envnames):
        envpaths.append(os.environ.get(envname))

    for i, envpath in enumerate(envpaths):
        if envpath is None:
            logmsg("error: environment variable {} needs to be declared".format(envpath))
            return -1

    elba_path = Path(envpaths[0]).resolve()
    combblas_path = Path(envpaths[1]).joinpath("CombBLAS").resolve()

    if not combblas_path.is_dir():
        logmsg("error: CombBLAS not found at '{}'".format(str(combblas_path)))
        return -1


    if argc < 2: return usage()

    compilation_path = Path("compilation").resolve()
    lower_kmer_bounds = [20]
    upper_kmer_bounds = [30]
    delta_chernoff = 0.1
    job_count = os.cpu_count()

    if osplat == "Linux":
        job_count //= 4

    try: opts, args = getopt.gnu_getopt(argv[1:], "L:U:d:I:h")
    except getopt.GetoptError as err:
        logmsg("error: {}".format(err))
        return usage()

    for o, a in opts:
        if o == "-h": return usage_extra()
        elif o == "-L": lower_kmer_bounds = list(map(lambda x: int(x), a.split(",")))
        elif o == "-U": upper_kmer_bounds = list(map(lambda x: int(x), a.split(",")))
        elif o == "-d": delta_chernoff = float(a)
        elif o == "-I": compilation_path = Path(a).resolve()

    if len(lower_kmer_bounds) != len(upper_kmer_bounds):
        logmsg("error: lower_kmer_bounds list '{}' and upper_kmer_bounds list '{}' are not the same length".format(str(lower_kmer_bounds), str(upper_kmer_bounds)))
        return usage()

    if len(args) != 1:
        return usage()

    progdir_path = Path(args[0]).resolve()
    progdir_path.mkdir(exist_ok=True)

    elba_build_path = compilation_path.joinpath("elba_build")
    elba_install_path = compilation_path.joinpath("elba_install")
    combblas_build_path = compilation_path.joinpath("combblas_build")
    combblas_install_path = compilation_path.joinpath("combblas_install")

    elba_build_path.mkdir(parents=True, exist_ok=True)
    elba_install_path.mkdir(parents=True, exist_ok=True)
    combblas_build_path.mkdir(parents=True, exist_ok=True)
    combblas_install_path.mkdir(parents=True, exist_ok=True)


    combblas_build_setup_cmd = ["cmake", "-S", str(combblas_path), "-B", str(combblas_build_path), "-DCMAKE_INSTALL_PREFIX=" + str(combblas_install_path)]
    if osplat == "Darwin": combblas_build_setup_cmd += ["-DCMAKE_C_COMPILER=gcc-11", "-DCMAKE_CXX_COMPILER=g++-11"]
    run_command(combblas_build_setup_cmd)

    combblas_build_cmd = ["make", "-j", str(job_count), "-C", str(combblas_build_path)]
    run_command(combblas_build_cmd)

    combblas_install_cmd = ["make", "install", "-j", str(job_count), "-C", str(combblas_build_path)]
    run_command(combblas_install_cmd)

    for l, u in zip(lower_kmer_bounds, upper_kmer_bounds):
        elba_build_setup_cmd = ["cmake", "-S", str(elba_path), "-B", str(elba_build_path), "-DLOWER_KMER_FREQ="+str(l), "-DUPPER_KMER_FREQ="+str(u), "-DDELTACHERNOFF="+str(delta_chernoff)]
        if osplat == "Darwin": elba_build_setup_cmd += ["-DCMAKE_C_COMPILER=gcc-11", "-DCMAKE_CXX_COMPILER=g++-11"]
        elba_build_cmd = ["make", "-j", str(job_count), "-C", str(elba_build_path)]

        run_command(elba_build_setup_cmd)
        run_command(elba_build_cmd)

        shutil.copy(str(elba_build_path.joinpath("elba")), str(progdir_path.joinpath("elba.l{}u{}".format(l, u))))

    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
