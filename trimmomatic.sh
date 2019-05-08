#!/bin/bash

#$ -N Trim
#$ -t 1-8
#$ -pe openmp 1

SAMPLE=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)

java -jar trimmomatic-0.36.jar \
PE -phred33 \
${SAMPLE}_1.fastq.gz \
${SAMPLE}_2.fastq.gz \
${SAMPLE}_1_paired.fastq.gz \
${SAMPLE}_1_unpaired.fastq.gz \
${SAMPLE}_2_paired.fastq.gz \
${SAMPLE}_2_unpaired.fastq.gz \
ILLUMINACLIP:data/users/aekimura/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36
