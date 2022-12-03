import sys
import getopt
import subprocess as sp
from pathlib import Path
import shutil
import os

def logmsg(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()

def usage():
    sys.stderr.write("Usage: python {} [options] /path/to/progdir\n\n".format(sys.argv[0]))
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -L INT,INT,...   prune k-mers below this bound [20]\n")
    sys.stderr.write("    -U INT,INT,...   prune k-mers above this bound [30]\n")
    sys.stderr.write("    -d FLOAT         chernoff bound [0.1]\n")
    sys.stderr.write("    -p PATH          ELBA main directory [../]\n")
    sys.stderr.write("    -I PATH          compilation directory [compilation]\n")
    return -1

def main(argc, argv):

    if argc < 2: return usage()

    source_path = Path("..").resolve()
    compilation_path = Path("compilation").resolve()
    lower_kmer_bounds = [20]
    upper_kmer_bounds = [30]
    delta_chernoff = 0.1

    try: opts, args = getopt.gnu_getopt(argv[1:], "L:U:d:p:I:h")
    except getopt.GetoptError as err:
        msglog("error: {}".format(err))
        return usage()

    for o, a in opts:
        if o == "-h": return usage()
        elif o == "-L": lower_kmer_bounds = list(map(lambda x: int(x), a.split(",")))
        elif o == "-U": upper_kmer_bounds = list(map(lambda x: int(x), a.split(",")))
        elif o == "-d": delta_chernoff = float(a)
        elif o == "-p": source_path = Path(a).resolve()
        elif o == "-I": compilation_path = Path(a).resolve()

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

    combblas_path = source_path.joinpath("CombBLAS")

    if not combblas_path.is_dir():
        logmsg("error: CombBLAS not found at '{}'".format(str(combblas_path)))

    combblas_build_setup_cmd = ["cmake", "-S", str(combblas_path), "-B", str(combblas_build_path), "-DCMAKE_INSTALL_PREFIX=" + str(combblas_install_path), "-DCMAKE_C_COMPILER=gcc-11", "-DCMAKE_CXX_COMPILER=g++-11"]
    logmsg(" ".join(combblas_build_setup_cmd))
    proc = sp.Popen(combblas_build_setup_cmd)
    proc.wait()

    combblas_build_cmd = ["make", "-j", str(12), "-C", str(combblas_build_path)]
    logmsg(" ".join(combblas_build_cmd))
    proc = sp.Popen(combblas_build_cmd)
    proc.wait()

    combblas_install_cmd = ["make", "install", "-j", str(12), "-C", str(combblas_install_path)]
    logmsg(" ".join(combblas_install_cmd))
    proc = sp.Popen(combblas_build_cmd)
    proc.wait()

    for l in lower_kmer_bounds:
        for u in upper_kmer_bounds:
            elba_build_setup_cmd = ["cmake", "-S", str(source_path), "-B", str(elba_build_path),
                                    "-DCMAKE_C_COMPILER=gcc-11", "-DCMAKE_CXX_COMPILER=g++-11",
                                    "-DLOWER_KMER_FREQ="+str(l), "-DUPPER_KMER_FREQ="+str(u),
                                    "-DDELTACHERNOFF="+str(delta_chernoff)]
            logmsg(" ".join(elba_build_setup_cmd))
            proc = sp.Popen(elba_build_setup_cmd)
            proc.wait()

            elba_build_cmd = ["make", "-j", str(12), "-C", str(elba_build_path)]
            logmsg(" ".join(elba_build_cmd))
            proc = sp.Popen(elba_build_cmd)
            proc.wait()

            shutil.copy(str(elba_build_path.joinpath("elba")), str(progdir_path.joinpath("elba.l{}u{}".format(l, u))))

    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
