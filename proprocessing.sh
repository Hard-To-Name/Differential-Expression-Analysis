#!/bin/bash

cat S288C_reference_sequence_R64-2-1_20150113.fsa | sed 's/>ref|NC_001133| [org=Saccharomyces cerevisiae] [str$

module load cufflinks
gffread saccharomyces_cerevisiae_R64-2-1_20150113.gff -T -o S288C.gtf

module load hisat2
hisat2-build -f S288C_reference_sequence_R64-2-1_20150113_new.fsa S288C
