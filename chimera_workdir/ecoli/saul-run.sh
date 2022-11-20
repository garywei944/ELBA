#!/bin/bash

export OMP_NUM_THREADS=1
# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads

NUM_READS=`grep '>' reads.fa | wc -l`

mpirun -np 1 ../../build_release/elba -i reads.fa -o ctg.fa --idxmap idmap -c $NUM_READS -k 31 --xa 15 -s 1 -O 100000 --afreq 100000 --mi -1 -g -1 -e -1 --ma 1 2>&1 | tee elba.out

minimap2 -x ava-pb -t12 reads.fa reads.fa > mm.paf

python ../../script/idname_replacement.py reads.fa elba.overlap.paf | sponge elba.overlap.paf
python ../../script/idname_replacement.py reads.fa elba.string.paf | sponge elba.string.paf
python ../../script/paf2gml.py reads.fa elba.overlap.paf elba.overlap.gml
python ../../script/paf2gml.py reads.fa elba.string.paf elba.string.gml
python ../../script/paf2gml.py reads.fa mm.paf mm.gml
