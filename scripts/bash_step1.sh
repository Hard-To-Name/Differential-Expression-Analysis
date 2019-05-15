#!/bin/bash

#$ -N Steps1-3
#$ -t 1-8
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

SAMPLE=$(head -n ${SGE_TASK_ID} ../data/samples.txt | tail -n 1)
GZ_DIR="../data/trimmomatic/"
GZ_1="_1_paired.fastq.gz"
GZ_2="_2_paired.fastq.gz"
REFERENCE=../S288C_reference_genome_R64-2-1_20150113/

if [ ! -d ../results ]; then
  mkdir -p ../results;

module load hisat2
hisat2 -p 1 --dta -x ${REFERENCE}index/S288C -1 ${GZ_DIR}${SAMPLE}${GZ_1} -2 ${GZ_DIR}${SAMPLE}${GZ_2} -S ../results/${SAMPLE}_S288C.sam


module load samtools
samtools sort -@ 1 -o ../results/${SAMPLE}_S288C.bam ../results/${SAMPLE}_S288C.sam


module load stringtie
stringtie -p 1 -G ${REFERENCE}/S288C.gtf -o ../results/${SAMPLE}_S288C.gtf -l ../results/${SAMPLE} ../results/${SAMPLE}_S288C.bam
