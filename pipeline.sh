#!/bin/bash
#
#$ -N Diff_Job
#$ -t 1-$NUM_OF_FILES
#$ -q free64,pub64 ### TBD
### -m beas
### -M duanr1@uci.edu
### -o
### -e
#$ -pe openmp 64-128  ##how many parallel environments to run and in what env i.e. openmp or mpi

JOB="/home/data/users/duanr1/bin/"
REFERENCE="/home/data/users/duanr1/bin/yeast"

###SEEDDIR="/home/data/users/duanr1/bin/RNAseq_data/"
###OUTDIR="/home/data/users/duanr1/bin/outputs/bt2/"
###IDXDIR="/home/data/users/duanr1/bin/drosophila_ref2/"

SEED=$(cat $JOBFILE | head -n $SGE_TASK_ID | tail -n 1)

###BOWTIE_OPTS=" -S "

###BTRUN="bowtie2 -x "${IDXDIR}"d_melanogaster_ffb5_22 "${SEEDDIR}${SEED}" -S "${OUTDIR}$(basename ${SEED} .txt)".sam"
$BTRUN

for FILE in *.fastq.gz; do

	hisat2 -p 8 --dta -x $REFERENCE/RNA_seq -1
	$REFERENCE/$FILE -2
	$REFERENCE/$FILE -S
	$JOB/${FILE%.fastq.gz}${SEED}.sam

done

for FILE in *.sam; do

	SAM = samtools sort -@ 8 -o "${FILE%.sam}.bam" $FILE

done

for FILE in *.bam; do

	stringtie -p 8 -G $REFERENCE/reference.gtf -o
	"${FILE%.bam}.gtf" -l $FILE

done

stringtie --merge -p 8 -G $REFERENCE/reference.gtf -o
stringtie_merged.gtf $REFERENCE/mergelist.txt

gffcompare -r $REFERENCE/reference.gtf -G -o merged stringtie_merge.gtf

stringtie ¨Ce ¨CB -p 8 -G stringtie_merged.gtf -o

for FILE in *.bam; do

	ballgown/${FILE%.bam}.gtf FILE

done