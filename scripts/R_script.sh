#!/bin/bash

#$ -N RScript
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

module load R/3.2.2

Rscript R_script_pipeline.R > ../results/R_output.txt
