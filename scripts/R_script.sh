#!/bin/bash

#$ -N RScript
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

source ~/.miniconda3testrc
conda activate r_env

Rscript rscript1.R
Rscript rscript2.R

conda deactivate
