#!/bin/bash

#$ -N Step6
#$ -t 1-8
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

REFERENCE=../S288C_reference_genome_R64-2-1_20150113/
SAMPLE=$(head -n ${SGE_TASK_ID} ../data/samples.txt | tail -n 1)

if [ ! -d ../ballgown ]; then
  mkdir -p ../ballgown;
fi

#module load stringtie
source ~/.miniconda3testrc
conda activate hisat2

stringtie -e -B -p 1 -G ${REFERENCE}stringtie_merged.gtf -o ../ballgown/${SAMPLE}/${SAMPLE}_S288C.gtf ../results/${SAMPLE}_S288C.bam

conda deactivate
