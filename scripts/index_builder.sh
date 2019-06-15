#!/bin/bash

#$ -N Index
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

FASTA_LOC=../S288C_reference_genome_R64-2-1_20150113/

# module load hisat2

source ~/.miniconda3testrc
conda activate hisat2

if [ ! -d ${FASTA_LOC}index ]; then
  mkdir -p ${FASTA_LOC}index;
fi

hisat2-build -f ${FASTA_LOC}S288C_reference_sequence_R64-2-1_20150113_new.fsa ${FASTA_LOC}index/S288C

conda deactivate
