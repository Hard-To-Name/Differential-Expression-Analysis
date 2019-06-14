#!/bin/sh

# This shell script generates an HTSeq file from (1) a raw FASTQ file containing the sequencing reads, (2) a FASTQ file containing 
# the corresponding index reads, and (3) a barcode text file.
# The loop in the main script at the end can be used to process many samples sequentially.
# Place this script into the folder with the original FASTQ files (before they have been multiplexed).
# Run this script by typing at the command prompt: ./fastq_pipeline.sh

# For this script, you need to install:
# Barcode Splitter: https://toolshed.g2.bx.psu.edu/repository?repository_id=7119c4f7a89efa57&changeset_revision=e7b7cdc1834d
# Samtools, which includes the gzip command: http://www.htslib.org/download/
# Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
# Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
# TopHat2: https://ccb.jhu.edu/software/tophat/index.shtml
# HTSeq: http://www-huber.embl.de/HTSeq/doc/overview.html

# Files describing the reference sequence and corresponding features of the genome are required.  For yeast, files are found here:
# http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/
# Download the compressed file "S288C_reference_genome_R64-2-1_20150113.tgz" and de-compress.
# Move "S288C_reference_sequence_R64-2-1_20150113.fsa" and "saccharomyces_cerevisiae_R64-2-1_20150113.gff" to a folder on your
# computer (e.g., "Yeast_genome"). 

# Bowtie2 needs to create index sequence files from FASTA file(s) containing the reference genome sequence.
# The following Bowtie2 function will generate these index files.
# Create the Yeast_genome directory before running this script.
cd /Users/hickmanm/Documents/Hypoxia_RNAseq/Yeast_genome
bowtie2-build -f S288C_reference_sequence_R64-2-1_20150113.fsa Scer_R64


# Download the compressed FASTQ files and the barcode file onto your computer into a folder that you create (e.g., FASTQ).
# Change to directory containing these files. De-compress these files with gzip.
# Read 1 file contains the actual sequence reads.  Read 2 file contains the Index (barcode) reads.
cd /Users/hickmanm/Documents/Hypoxia_RNAseq/FASTQ
gzip -d read1.fastqsanger.gz
gzip -d read2_index.fastqsanger.gz


# A Barcode splitter is used to de-multiplex the original FASTQ file (i.e., generated from multiplexing several samples into one sequencing lane),
# producing one FASTQ file for each sample. The resulting FASTQ files are submitted to GEO, and also used for subsequent steps in this script.
# Note that sometimes, sequencing will result in a pre-multiplex sequencing FASTQ file (read 1) and a separate index FASTQ file (read 2);
# this is in contrast to a single FASTQ file containing both the sequencing reads and corresponding index reads.  The version of the Barcode
# splitter on Galaxy does not allow two separate FASTQ files.  Lance Parsons at Princeton University created a PERL version of Barcode splitter
# which is used here that inputs two FASTQ files. This version can be found at:
# https://toolshed.g2.bx.psu.edu/repository?repository_id=7119c4f7a89efa57&changeset_revision=e7b7cdc1834d
# The file barcode.txt is a tab-delimited file that pairs a nucleotide barcode with its corresponding sample.
cat read1.fastqsanger | perl /usr/local/bin/fastx_barcode_splitter.pl --bcfile barcodes.txt --prefix /Users/hickmanm/Documents/Hypoxia_RNAseq/FASTQ/ --suffix ".fastq"  --idxfile read2_index.fastqsanger --idxidstrip 1  --mismatches 0



## Define the get_count function to be carried out on each FASTQ file that has been created by de-multiplexing.
# Note that the $1 variable holds the input (e.g., NH_42_1 during the first loop). Also, this function assumes certain directories
# which should be set up on your computer. You may have to change the file paths below, depending on your directory names.
get_count() {

# Change to directory containing the FASTQ files which were just made from de-multiplexing.
# Create this directory before running this script.
cd /Users/hickmanm/Documents/Hypoxia_RNAseq/FASTQ

# Quality Trim the reads from the 5' and 3' ends using Trimmomatic with default settings.
# Note: fastq_quality_trimmer from the FASTX toolkit is not sufficient, because it only trims one end of the read.
java -jar /Applications/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 $1.fastq $1_trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Map the reads to the reference genome using TopHat, which can map across introns. TopHat uses the Bowtie2 algorithm to help
# with mapping. (Note that BWA and Bowtie, working alone, do not allow mapping across introns.)
tophat --library-type fr-firststrand  /Users/hickmanm/Documents/Hypoxia_RNAseq/Yeast_genome/Scer_R64 $1_trim.fastq

# Move the TopHat output into a TopHat directory. Create this directory before running this script. 
cd tophat_out
mv accepted_hits.bam /Users/hickmanm/Documents/Hypoxia_RNAseq/TopHat/$1.bam
cd ..
rm -R tophat_out

# Change to directory containing the BAM files created by TopHat. Create this directory before running this script.
cd /Users/hickmanm/Documents/Hypoxia_RNAseq/TopHat

# Count the number of reads mapping to each S. cerevisiae gene, using an HTSeq command. The gff file contains the locations of the genome features
# (e.g., genes). The output of this command is piped into a tab-delimited text file. This file can be submitted to GEO as raw read counts for the sample.
htseq-count -f bam -s yes -t gene -i ID $1.bam /Users/hickmanm/Documents/Hypoxia_RNAseq/Yeast_genome/saccharomyces_cerevisiae_R64-2-1_20150113.gff > /Users/hickmanm/Documents/Hypoxia_RNAseq/HTSEQ/$1.txt

}




## MAIN SCRIPT ##

# This script consists of a loop to be repeated for each base string provided (e.g., NH_42_1). For example, in the first loop, NH_42_1 will be input
# for the get_count function defined above. In this function, NH_42_1.fastq will be trimmed to generate NH_42_1_trim.fastq. This trimmed file
# will be mapped using Tophat to generate NH_42_1.bam. This BAM file will be input for htseq-count, thus generating NH_42_1.txt.
# Thus, you should name your original FASTQ files accordingly.
for file in NH_42_1 NH_42_2 NH_42_3 NH_42_4 NH_42_5 NH_42_6 NH_42_7 NH_42_8
do
 get_count $file
done