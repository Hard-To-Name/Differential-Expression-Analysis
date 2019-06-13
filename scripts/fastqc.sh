#!/bin/bash

#$ -N FastQC
#$ -t 1-8
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

FASTQ_LOC=../data/
SAMPLE=$(head -n ${SGE_TASK_ID} ../data/samples.txt | tail -n 1)

conda activate hisat2

fastqc ${FASTQ_LOC}${SAMPLE}_1.fastq.gz
fastqc ${FASTQ_LOC}${SAMPLE}_2.fastq.gz

if [ ! -d ${FASTQ_LOC}fastqc ]; then
  mkdir -p ${FASTQ_LOC}fastqc;
fi

mv *.html ${FASTQ_LOC}fastqc
mv *.zip ${FASTQ_LOC}fastqc

conda deactivate