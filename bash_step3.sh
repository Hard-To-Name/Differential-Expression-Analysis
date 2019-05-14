#!/bin/bash

#$ -N Step6
#$ -t 1-8
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

SAMPLE=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)

module load stringtie
stringtie -e -B -p 1 -G stringtie_merged.gtf -o ballgown/${SAMPLE}/${SAMPLE}_S288C.gtf ${SAMPLE}_S288C.bam
