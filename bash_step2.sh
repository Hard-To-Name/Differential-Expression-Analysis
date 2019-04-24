#!/bin/bash
#$ -N Step_1
#$ -t 1-12
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1
#$ -R y

SEED=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)
GZ_DIR="chrX_data/samples/"
GTF="chrX_data/genes/chrX.gtf"
GZ_1="_chrX_1.fastq.gz"
GZ_2="_chrX_2.fastq.gz"
GFF_COMPARE="/data/users/duanr1/bin/gffcompare/gffcompare"
BGN_DATA_DIR="ballgown/"

module load stringtie

stringtie --merge -p 1 -G ${GTF} -o stringtie_merged.gtf chrX_data/mergelist.txt

${GFF_COMPARE} -r ${GTF} -G -o merged stringtie_merged.gtf

stringtie -e -B -p 1 -G stringtie_merged.gtf -o ${BGN_DATA_DIR}${SEED}/${SEED}_chrX.gtf ${SEED}_chrX.bam
~
