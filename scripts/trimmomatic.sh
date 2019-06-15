#!/bin/bash

#$ -N Trim
#$ -t 1-8
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

source ~/.miniconda3testrc
conda activate trimmomatic

FASTQ_LOC=../data/
SAMPLE=$(head -n ${SGE_TASK_ID} ${FASTQ_LOC}samples.txt | tail -n 1)

if [ ! -d ${FASTQ_LOC}trimmomatic ]; then
  mkdir -p ${FASTQ_LOC}trimmomatic;
fi

java -jar ../Trimmomatic-0.36/trimmomatic-0.36.jar \
PE -phred33 \
${FASTQ_LOC}${SAMPLE}_1.fastq.gz \
${FASTQ_LOC}${SAMPLE}_2.fastq.gz \
${FASTQ_LOC}trimmomatic/${SAMPLE}_1_paired.fastq.gz \
${FASTQ_LOC}trimmomatic/${SAMPLE}_1_unpaired.fastq.gz \
${FASTQ_LOC}trimmomatic/${SAMPLE}_2_paired.fastq.gz \
${FASTQ_LOC}trimmomatic/${SAMPLE}_2_unpaired.fastq.gz \
ILLUMINACLIP:../Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:0 \
TRAILING:0 \
SLIDINGWINDOW:4:15 \
MINLEN:5 \
AVGQUAL:20

conda deactivate
