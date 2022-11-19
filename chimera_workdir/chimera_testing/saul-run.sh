#!/bin/bash

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

srun -n 1 -c 256 ../../build_release/elba -i reads.fa -o ctg.fa --idxmap idmap -c 206 -k 31 --xa 15 -s 1 -O 100000 --afreq 100000 --ch 0 --mi -1 -g -1 -e -1 --ma 1 2>&1 | tee elba.out

