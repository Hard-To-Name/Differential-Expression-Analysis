#!/bin/bash

#$ -N Steps4-5
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

REFERENCE="../S288C_reference_genome_R64-2-1_20150113/"

#module load stringtie

source ~/.miniconda3testrc
conda activate hisat2

stringtie --merge -p 1 -G ${REFERENCE}S288C.gtf -o ${REFERENCE}stringtie_merged.gtf ${REFERENCE}mergelist.txt

gffread -E ${REFERENCE}stringtie_merged.gtf -o- > ${REFERENCE}stringtie_merged.gff3 # for bedtools use

/data/users/duanr1/bin/gffcompare/gffcompare -r ${REFERENCE}S288C.gtf -G -o ${REFERENCE}merged ${REFERENCE}stringtie_merged.gtf

conda deactive
