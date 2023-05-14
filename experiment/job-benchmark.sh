#!/bin/bash
#SBATCH -N 1
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

# MPI processes scale
for n in 25 16 9 4 1; do
    # injection
    srun -N 1 -n "$n" -c 4 --cpu_bind=cores -G 4 --gpu-bind=none \
        ../build_release/elba -i \
        "$DATA" \
        -k 31 --idxmap ecoli-idxmap -c "$SEQUENCE_COUNT" --alph dna \
        --af ecoli-gpu -s 1 -O 100000 --afreq 100000 --ga 15 --gpu_num 4 |&
        tee outputs/"$NAME"/injection_N_1_n_"$n"_gpu_4.out

    # one-to-all
    srun -N 1 -n "$n" -c 4 --cpu_bind=cores -G 4 --gpu-bind=none \
        ../build_release/elba -i \
        "$DATA" \
        -k 31 --idxmap ecoli-idxmap -c "$SEQUENCE_COUNT" --alph dna \
        --af ecoli-gpu -s 1 -O 100000 --afreq 100000 --ga 15 --gpu_num -4 |&
        tee outputs/"$NAME"/one2all_N_1_n_"$n"_gpu_4.out
done

for gpu in {1..2}; do
    # injection
    srun -N 1 -n 16 -c 4 --cpu_bind=cores -G "$gpu" --gpu-bind=none \
        ../build_release/elba -i \
        "$DATA" \
        -k 31 --idxmap ecoli-idxmap -c "$SEQUENCE_COUNT" --alph dna \
        --af ecoli-gpu -s 1 -O 100000 --afreq 100000 --ga 15 \
        --gpu_num "$gpu" |&
        tee outputs/"$NAME"/injection_N_1_n_16_gpu_"$gpu".out

    # one-to-all
    srun -N 1 -n 16 -c 4 --cpu_bind=cores -G "$gpu" --gpu-bind=none \
        ../build_release/elba -i \
        "$DATA" \
        -k 31 --idxmap ecoli-idxmap -c "$SEQUENCE_COUNT" --alph dna \
        --af ecoli-gpu -s 1 -O 100000 --afreq 100000 --ga 15 \
        --gpu_num -"$gpu" |&
        tee outputs/"$NAME"/one2all_N_1_n_16_gpu_"$gpu".out
done
