#!/bin/bash


mkdir S288C_data

cd S288C_data

mkdir genes
mkdir genome
mkdir indexes
mkdir samples

cd ../

mv *_new.fsa /S288C_data/genome
mv *.ht2 /S288C_data/indexes
mv *_paired.fastq.gz /S288C_data/samples
mv S228C.gtf /S288C_data/genes
mv mergelist.txt S288C_data
mv S228C_phenodata.csv S288C_data
