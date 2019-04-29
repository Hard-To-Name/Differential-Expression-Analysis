#!/bin/bash
#$ -N Step_1
#$ -t 1-12
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1
#$ -R y

SEED=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)
BGN_DATA_DIR="ballgown/"

module load stringtie

stringtie -e -B -p 1 -G stringtie_merged.gtf -o ${BGN_DATA_DIR}${SEED}/${SEED}_chrX.gtf ${SEED}_chrX.bam
