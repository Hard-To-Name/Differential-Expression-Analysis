#!/bin/bash
#$ -N Step_1
#$ -ckpt restart
#$ -q bio,abio*,pub64,free64,epyc
#$ -pe openmp 1

if [ ! -d ../data ]; then
  mkdir -p ../data;
fi

DATA_LOC=../data/

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/000/SRR1257640/SRR1257640_1.fastq.gz -o ${DATA_LOC}SRR1257640_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/000/SRR1257640/SRR1257640_2.fastq.gz -o ${DATA_LOC}SRR1257640_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/003/SRR1257793/SRR1257793_1.fastq.gz -o ${DATA_LOC}SRR1257793_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/003/SRR1257793/SRR1257793_2.fastq.gz -o ${DATA_LOC}SRR1257793_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/007/SRR1259267/SRR1259267_1.fastq.gz -o ${DATA_LOC}SRR1259267_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/007/SRR1259267/SRR1259267_2.fastq.gz -o ${DATA_LOC}SRR1259267_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/005/SRR1514795/SRR1514795_1.fastq.gz -o ${DATA_LOC}SRR1514795_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/005/SRR1514795/SRR1514795_2.fastq.gz -o ${DATA_LOC}SRR1514795_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/005/SRR1515155/SRR1515155_1.fastq.gz -o ${DATA_LOC}SRR1515155_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/005/SRR1515155/SRR1515155_2.fastq.gz -o ${DATA_LOC}SRR1515155_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/007/SRR1257637/SRR1257637_1.fastq.gz -o ${DATA_LOC}SRR1257637_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/007/SRR1257637/SRR1257637_2.fastq.gz -o ${DATA_LOC}SRR1257637_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/006/SRR1515156/SRR1515156_1.fastq.gz -o ${DATA_LOC}SRR1515156_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/006/SRR1515156/SRR1515156_2.fastq.gz -o ${DATA_LOC}SRR1515156_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/004/SRR1514794/SRR1514794_1.fastq.gz -o ${DATA_LOC}SRR1514794_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/004/SRR1514794/SRR1514794_2.fastq.gz -o ${DATA_LOC}SRR1514794_2.fastq.gz
