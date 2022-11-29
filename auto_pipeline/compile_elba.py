import sys
import getopt
import subprocess as sp
from pathlib import Path
import shutil
import os

def usage():
    sys.stderr.write("Usage: python {} [options] <elba>\n\n".format(sys.argv[0]))
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -L INT   prune k-mers below this bound [20]\n")
    sys.stderr.write("    -U INT   prune k-mers above this bound [30]\n")
    sys.stderr.write("    -d FLOAT chernoff bound [0.1]\n")
    sys.stderr.write("    -p FILE  relative path to ELBA main directory [../]\n")
    sys.stderr.write("    -M       on personal mac computer\n")
    return -1

def main(argc, argv):

    if argc < 2: return usage()

    elba_pathname = "../"
    lower_kmer_bound = 20
    upper_kmer_bound = 30
    delta_chernoff = 0.1
    on_mac=False

    try: opts, args = getopt.gnu_getopt(argv[1:], "L:U:d:p:Mh")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        return usage()

    for o, a in opts:
        if o == "-h": return usage()
        elif o == "-L": lower_kmer_bound = int(a)
        elif o == "-U": upper_kmer_bound = int(a)
        elif o == "-d": delta_chernoff = float(a)
        elif o == "-p": elba_path = a
        elif o == "-M": on_mac = True

    if len(args) != 1:
        return usage()

    exe_path = Path(args[0])

    elba_path = Path(elba_pathname).resolve()
    build_path = elba_path.joinpath("build_release")
    combblas_path = elba_path.joinpath("CombBLAS").resolve()

    build_path.mkdir(exist_ok=True)
    build_path = build_path.resolve()

    combblas_build_path = combblas_path.joinpath("build").resolve()
    combblas_install_path = combblas_path.joinpath("install").resolve()

    combblas_build_path.mkdir(exist_ok=True)
    combblas_install_path.mkdir(exist_ok=True)


    combblas_cmake_cmd = ["cmake", "-S", str(combblas_path.absolute()), "-B", str(combblas_build_path.absolute()), "-DCMAKE_INSTALL_PREFIX={}".format(str(combblas_install_path.absolute()))]

    if on_mac: combblas_cmake_cmd += ["-DCMAKE_CXX_COMPILER=g++-11", "-DCMAKE_C_COMPILER=gcc-11"]

    sys.stdout.write(" ".join(combblas_cmake_cmd) + "\n")
    proc = sp.Popen(combblas_cmake_cmd)
    proc.wait()

    job_count = 64
    if on_mac: job_count = 12

    combblas_build_cmd = ["make", "-j", str(job_count), "-C", str(combblas_build_path.absolute())]
    sys.stdout.write(" ".join(combblas_build_cmd) + "\n")
    proc = sp.Popen(combblas_build_cmd)
    proc.wait()

    combblas_install_cmd = ["make", "install", "-j", str(job_count), "-C", str(combblas_build_path.absolute())]
    sys.stdout.write(" ".join(combblas_install_cmd) + "\n")
    proc = sp.Popen(combblas_install_cmd)
    proc.wait()

    elba_cmake_cmd = ["cmake", "-S", str(elba_path.absolute()), "-B", str(build_path.absolute()),
                      "-DLOWER_KMER_FREQ={}".format(lower_kmer_bound),
                      "-DUPPER_KMER_FREQ={}".format(upper_kmer_bound),
                      "-DDELTACHERNOFF={}".format(delta_chernoff)]

    if on_mac: elba_cmake_cmd += ["-DCMAKE_CXX_COMPILER=g++-11", "-DCMAKE_C_COMPILER=gcc-11"]

    sys.stdout.write(" ".join(elba_cmake_cmd) + "\n")
    proc = sp.Popen(elba_cmake_cmd)
    proc.wait()

    elba_build_cmd = ["make", "-j", str(job_count), "-C", str(build_path.absolute())]
    sys.stdout.write(" ".join(elba_build_cmd) + "\n")
    proc = sp.Popen(elba_build_cmd)
    proc.wait()

    shutil.copy(str(build_path.joinpath("elba").absolute()), str(exe_path.absolute()))

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
