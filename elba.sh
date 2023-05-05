#!/bin/bash

#source corigpu-env.sh
export COMBBLAS_HOME=$PWD
export BLOOM_HOME=$PWD/src/libbloom/  # k-mer counting phase
export SEQAN_HOME=$PWD/seqan    # cpu based aligner
export LOGAN_HOME=$PWD/LoganGPU # gpu based aligner


cd $COMBBLAS_HOME
#rm -rf build_release
#mkdir build_release # elba executable in here
cd build_release
cmake -DLOWER_KMER_FREQ=20 -DUPPER_KMER_FREQ=30 -DDELTACHERNOFF=0.1 ..
make -j8

echo ""
echo "ELBA installation completed."
echo ""