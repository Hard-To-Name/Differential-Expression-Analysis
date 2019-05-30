#!/bin/bash

#$ -N Gffread
#$ -t 1
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

FASTA_LOC=../S288C_reference_genome_R64-2-1_20150113/

source ~/.miniconda3testrc
conda activate hisat2 

gffread ${FASTA_LOC}saccharomyces_cerevisiae_R64-2-1_20150113.gff -T -o ${FASTA_LOC}S288C.gtf

conda deactivate
