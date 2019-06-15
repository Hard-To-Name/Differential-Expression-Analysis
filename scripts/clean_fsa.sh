#!/bin/bash

#$ -N CleanFSA
#$ -t 1
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1


cat ../S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa | sed 's/>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]/I/g' | awk '/^>/{print ">" ++i; next}{print}' > ../S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113_new.fsa

