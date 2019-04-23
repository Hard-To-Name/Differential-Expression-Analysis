#!/bin/bash
#$ -N Step_1
#$ -t 1-12
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1
#$ -R y

SEED=$(head -n ${SGE_TASK_ID} samples.txt | tail -n 1)
GZ_DIR="chrX_data/samples/"
GZ_1="_chrX_1.fastq.gz"
GZ_2="_chrX_2.fastq.gz"


###hisat2 -p 1 --dta -x chrX_data/indexes/chrX_tran -1 ${GZ_DIR}${SAMPLE}${GZ_1} -2 ${GZ_DIR}${SAMPLE}${GZ_2} -S ${SAMPL$

echo ${GZ_DIR}${SEED}${GZ_1}
