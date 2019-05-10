#!/bin/bash

#$ -N Index
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

module load hisat2

hisat2-build -f S288C_reference_sequence_R64-2-1_20150113.fasta S288C
