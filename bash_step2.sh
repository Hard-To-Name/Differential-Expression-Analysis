#!/bin/bash
#$ -N Step_1
#$ -t 1
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1
#$ -R y

GZ_DIR="chrX_data/samples/"
GTF="chrX_data/genes/chrX.gtf"
GFF_COMPARE="/data/users/duanr1/bin/gffcompare/gffcompare"

module load stringtie

stringtie --merge -p 1 -G ${GTF} -o stringtie_merged.gtf chrX_data/mergelist.txt

${GFF_COMPARE} -r ${GTF} -G -o merged stringtie_merged.gtf