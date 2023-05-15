#!/bin/bash
#SBATCH -N 2
#SBATCH -C gpu
#SBATCH -G 4
#SBATCH -q regular
#SBATCH -J ecoli-cpu-gpu
#SBATCH -t 02:00:00
#SBATCH -A m4341
#SBATCH --error=logs/ecoli-elba-%j
#SBATCH --output=logs/ecoli-elba-%j
#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=true

# salloc --nodes 1 --gpus 4 --qos interactive --time 02:00:00 --constraint gpu --account=m4341

# # 29x dataset
# DATA=../ecoli_hifi_29x.fasta
# NAME=ecoli29x
# SEQUENCE_COUNT=8605

# 100x dataset
DATA=../ecoli100x.fa
NAME=ecoli100x
SEQUENCE_COUNT=91394

mkdir -p outputs/$NAME

run() {
    gpu_num=$gpu

    srun -N "$N" --ntasks-per-node "$n" \
        -c 4 --cpu_bind=cores \
        --gpus-per-node "$gpu" --gpu-bind=none \
        ../build_release/elba -i "$DATA" \
        -k 31 --idxmap ecoli-idxmap -c "$SEQUENCE_COUNT" --alph dna \
        --af ecoli-gpu -s 1 -O 100000 --afreq 100000 --ga 15 \
        --gpu_num "$gpu_num" |&
        tee outputs/"$NAME"/"$job"_N_"$N"_n_"$n"_gpu_"$gpu".out
}

job=opt_injection

# N=1
# for n in 25 16 9 4 1; do
#     gpu=4
#     run
# done
# for gpu in 1 2; do
#     n=16
#     run
# done

N=2
for n in 2 8 18; do
    gpu=4
    run
done
for gpu in 1 2; do
    n=8
    run
done

rm *.txt
rm elba.contigs.fa
rm *.mtx
rm ecoli-gpu-readnamemap*
