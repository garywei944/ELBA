#!/bin/bash

export COMBBLAS_HOME=$PWD
export ELBA_HOME=$PWD
export BLOOM_HOME=$PWD/src/libbloom/
export SEQAN_HOME=$PWD/seqan
export OMP_NUM_THREADS=1
export PATH=$PATH:$ELBA_HOME/auto_pipeline_modules
