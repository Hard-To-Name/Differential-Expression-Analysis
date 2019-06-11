#!/bin/bash

#$ -N Bedtools
#$ -t 1-10
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

source ~/.miniconda3testrc
conda activate bedtools

REFERENCE=../S288C_reference_genome_R64-2-1_20150113/
SIGNIFICANT=$(head -n ${SGE_TASK_ID} ${REFERENCE}significant1.txt | tail -n 1)

echo ${SIGNIFICANT}.gtf

if [ ! -d ${REFERENCE}significant ]; then
  mkdir -p ${REFERENCE}significant;
fi

grep -a ${SIGNIFICANT} ${REFERENCE}stringtie_merged.gtf > ${REFERENCE}significant/${SIGNIFICANT}.gtf

bedtools getfasta -fi ${REFERENCE}S288C_reference_sequence_R64-2-1_20150113_new.fsa -bed ${REFERENCE}significant/${SIGNIFICANT}.gtf -fo ${REFERENCE}significant/${SIGNIFICANT}.fasta

sed -i '/>/d' ${REFERENCE}significant/${SIGNIFICANT}.fasta

conda deactivate
