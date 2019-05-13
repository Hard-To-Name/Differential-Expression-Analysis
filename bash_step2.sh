#!/bin/bash

#$ -N Steps4-5
#$ -t 1-8
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1


module load stringtie
stringtie --merge -p 1 -G S288C_data/genes/S288C.gtf -o stringtie_merged.gtf  S288C_data/mergelist.txt

modlue load gffcompare
gffcompare -r S288C_data/genes/S288C.gtf -G -o merged stringtie_merged.gtf
