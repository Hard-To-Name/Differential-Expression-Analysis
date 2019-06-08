#!/bin/bash

#$ -N Bedtools
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

source ~/.miniconda3testrc
conda activate bedtools

REFERENCE=../S288C_reference_genome_R64-2-1_20150113/

bedtools getfasta -fi ${REFERENCE}S288C_reference_sequence_R64-2-1_20150113.fsa -bed ${REFERENCE}stringtie_merged.gff3 -fo ${REFERENCE}bedtools_output.fasta

conda deactivate
