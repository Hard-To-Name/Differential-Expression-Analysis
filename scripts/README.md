# Background

## Global Variables

```sh
FASTQ_LOC=../data/
FASTA_LOC=../S288C_reference_genome_R64-2-1_20150113/
GZ_DIR="../data/trimmomatic/"
GZ_1="_1_paired.fastq.gz"
GZ_2="_2_paired.fastq.gz"
REFERENCE=../S288C_reference_genome_R64-2-1_20150113/
SAMPLE=$(head -n ${SGE_TASK_ID} ../data/samples.txt | tail -n 1)
SIGNIFICANT=$(head -n ${SGE_TASK_ID} ${REFERENCE}significant1.txt | tail -n 1)

```

## Directory Structure
```
.  
├── scripts  
|   ├── pipeline.sh # include every other bash scripts in order  
|   └── *.sh  
├── data  
|   ├── fastqc  
|   |   ├── *.html  
|   |   └── *.zip  
|   └── trimmomatic  
|       ├── *_paired.fastq.gz  
|       └── *_unpaired.fastq.gz  
├── S288C_reference_genome_R64-2-1_20150113  
|   ├── index  
|   |   └── *.ht2  
|   ├── significant1 (for 4 samples in group 30°C)  
|   |   ├── *.gtf  
|   |   └── *.fasta  
|   ├── significant2 (for 4 samples in group 37°C)  
|   |   ├── *.gtf  
|   |   └── *.fasta  
|   ├── .gtf of merged annotation from StringTie
|   ├── .gtf of reference genome  
|   ├── .fsa of reference genome before and after cleaning headers  
|   ├── phenotype.csv for ballgown  
|   └── merge information  
├── results  
|   ├── *.sam before merge  
|   ├── *.bam before merge  
|   └── *.gtf before merge  
└── ballgown  
    ├── temp30 (data for temperature 30)  
    └── temp37 (data for temperature 37)  
```

## Environment

### Conda Information
* Miniconda3 version 4.5.12
* Python version 3.7.1.final.0
* Platform: linux-64
  
Configuration File: .miniconda3testrc
```sh
# added by Miniconda3 4.5.12 installer
# >>> conda init >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$(CONDA_REPORT_ERRORS=false '/data/users/duanr1/bin/miniconda3/bin/conda' shell.bash hook 2> /dev/null)"
if [ $? -eq 0 ]; then
    \eval "$__conda_setup"
else
    if [ -f "/data/users/duanr1/bin/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data/users/duanr1/bin/miniconda3/etc/profile.d/conda.sh"
        CONDA_CHANGEPS1=false conda activate base
    else
        \export PATH="/data/users/duanr1/bin/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda init <<<
```

### Conda Environments
* hisat2
Configuration file [hisat2_env.yml](conda_config/hisat2_env.yml)

* trimmomatic
Configuration file [trimmomatic_env.yml](conda_config/trimmomatic_env.yml)

* r_env
Configuration file [r_env.yml](conda_config/r_env.yml)

* bedtools
Configuration file [bedtools_env.yml](conda_config/bedtools_env.yml)  
  
Instructions of installing Miniconda3 can be found [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Environments can be created by the following command:  
```conda env create -f environment.yml```

***

# Pipeline

## Downloading Data

### Downloading Sample Data
We are using Saccharomyces cerevisiae S228C sample reads from the European Bioinformatics Institute.  
Download all 8 samples by executing ```./sccripts/download_fastq.sh```

### Downloading Reference Data
We are using Saccharomyces cerevisiae S228C Release R64-2-1 from  
https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/  
Download ```S288C_reference_genome_R64-2-1_20150113.tgz``` and unzip it to ```./S288C_reference_genome_R64-2-1_20150113```

## Preprocessing

#### Quality Measurement of Sample Reads
Use FASTQC to measure the quality scores (phred) across all bases of each sample. Results are available at [here](../data/fastqc).
```sh
fastqc ${FASTQ_LOC}${SAMPLE}_1.fastq.gz
fastqc ${FASTQ_LOC}${SAMPLE}_2.fastq.gz

if [ ! -d ${FASTQ_LOC}fastqc ]; then
  mkdir -p ${FASTQ_LOC}fastqc;
fi

mv *.html ${FASTQ_LOC}fastqc
mv *.zip ${FASTQ_LOC}fastqc
```

#### Cleaning Reference Genome
Remove unnecessary headers of all .fastq files. The complete script is available [here](clean_fsa.sh).
```sh
cat ../S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa | sed 's/>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]/I/g' | awk '/^>/{print ">" ++i; next}{print}' > ../S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113_new.fsa
```

#### Converting gff3 to gtf
Under conda environment hisat2, convert reference genome annotation file in .gff3 format to ballgown required format .gtf. The complete script is available [here](gffread.sh).

```sh
gffread ${FASTA_LOC}saccharomyces_cerevisiae_R64-2-1_20150113.gff -T -o ${FASTA_LOC}S288C.gtf
```

#### Building Index
Under conda environment hisat2, build HISAT2 index for the reference genome. The complete script is available [here](index_builder.sh).

```sh
if [ ! -d ${FASTA_LOC}index ]; then
  mkdir -p ${FASTA_LOC}index;
fi

hisat2-build -f ${FASTA_LOC}S288C_reference_sequence_R64-2-1_20150113_new.fsa ${FASTA_LOC}index/S288C
```

#### Trimming Sample Reads
Under conda environment trimmomatic, generate both paired reads and unpaired reads at ```./data/trimmomatic``` (Only paired data are used). The complete script is available [here](trimmomatic.sh).

```sh
if [ ! -d ${FASTQ_LOC}trimmomatic ]; then
  mkdir -p ${FASTQ_LOC}trimmomatic;
fi

java -jar ../Trimmomatic-0.36/trimmomatic-0.36.jar \
PE -phred33 \
${FASTQ_LOC}${SAMPLE}_1.fastq.gz \
${FASTQ_LOC}${SAMPLE}_2.fastq.gz \
${FASTQ_LOC}trimmomatic/${SAMPLE}_1_paired.fastq.gz \
${FASTQ_LOC}trimmomatic/${SAMPLE}_1_unpaired.fastq.gz \
${FASTQ_LOC}trimmomatic/${SAMPLE}_2_paired.fastq.gz \
${FASTQ_LOC}trimmomatic/${SAMPLE}_2_unpaired.fastq.gz \
ILLUMINACLIP:../Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:0 \
TRAILING:0 \
SLIDINGWINDOW:4:15 \
MINLEN:5 \
AVGQUAL:20
```

## Results

***

# Pipeline Steps

## Description

**1. Map the reads for each sample to the reference genome**  
Use hisat2 to:  
(1) build indices of sample reads  
(2) map sample reads to the reference genome  
(3) output .sam files  

```sh
hisat2 -p 1 --dta -x ${REFERENCE}index/S288C -1 ${GZ_DIR}${SAMPLE}${GZ_1} -2 ${GZ_DIR}${SAMPLE}${GZ_2} -S ../results/${SAMPLE}_S288C.sam
```
  
**2. Sort and convert the SAM files to BAM**  
Use samtools to sort .sam files to .bam files such that:  
(1) the files can store alignment information generated by various alignment programs  
(2) the files are simple to be generated and converted  
(3) the files are compact in file size  
(4) the files are indexed by genomic position to efficiently retrieve all reads aligning to a locus  

```sh
samtools sort -@ 1 -o ../results/${SAMPLE}_S288C.bam ../results/${SAMPLE}_S288C.sam
```
  
**3. Assemble transcripts for each sample**  
Use stringtie to assemble RNA-Seq alignments into potential transcripts and output an estimated .gtf file for each sample  

```sh
stringtie -p 1 -G ${REFERENCE}S288C.gtf -o ../results/${SAMPLE}_S288C.gtf -l ../results/${SAMPLE} ../results/${SAMPLE}_S288C.bam
```
  
**4. Merge transcripts from all samples**  
Use stringtie to merge all transcripts generated by the previous step into a non-redundant set of transcripts and output a merged .gtf file  

```sh
stringtie --merge -p 1 -G ${REFERENCE}S288C.gtf -o ${REFERENCE}stringtie_merged.gtf ${REFERENCE}mergelist.txt
```

**5. Examine how the transcripts compare with the reference annotation (optional)**  
Use gffcompare to compare and compute the estimated accuracy of the merged .gtf file and .gtf files of each sample, to check how the predicted transcripts relate to an annotation file  

```sh
gffcompare -r ${REFERENCE}S288C.gtf -G -o ${REFERENCE}merged ${REFERENCE}stringtie_merged.gtf
```
  
**6. Estimate transcript abundances and create table counts for Ballgown**  
Use stringtie to generate .gtf file and .bam file for each sample reads based on the merged annotation file

```sh
stringtie -e -B -p 1 -G ${REFERENCE}stringtie_merged.gtf -o ../ballgown/${SAMPLE}/${SAMPLE}_S288C.gtf ../results/${SAMPLE}_S288C.bam
```

There will be 8 files for the samples in ```../ballgown```. Create two files named ```temp30``` and ```temp37``` respectively; then put all 4 files listed in ```phenotype1.csv``` in folder ```temp30``` and put all 4 files listed in ```phenotype2.csv``` in folder ```temp37```.

**7. Perform statistical analysis**
The script for group under temperature 30°C is available [here](rscript1.R); The script for group under temperature 37°C is available [here](rscript2.R);

* Load relevant R packages
```R
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
```

* Load the phenotype data for the samples
```R
pheno_data = read.csv("../S288C_reference_genome_R64-2-1_20150113/phenotype1.csv")
```

* Read in te expression data that were calculated by StringTie.
```R
bg_genome = ballgown(dataDir = "../ballgown/temp30", samplePattern = "SRR", pData = pheno_data)
```

* Filter to remove low-abundance genes
```R
bg_genome_filt = subset(bg_genome, "rowVars(texpr(bg_genome)) > 1", genomesubset=TRUE)
```

* Identify transcripts that show statistically significant differences between groups
```R
results_transcripts = stattest(bg_genome_filt, feature="transcript",covariate="population", getFC=TRUE, meas="FPKM")
```

* Identify genes that show statistically significant differences between groups
```R
results_genes = stattest(bg_genome_filt, feature="gene", covariate="population", getFC=TRUE, meas="FPKM")
```

* Add gene names and gene IDs to the results_transcripts data frame
```R
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_genome_filt), geneIDs=ballgown::geneIDs(bg_genome_filt), results_transcripts)
```

* Sort the results from the smallest P value to the largest; Select top 10 transcripts with the smallest p values; Write the results to ```../S288C_reference_genome_R64-2-1_20150113/``` in a .csv format.
```R
results_transcripts = arrange(results_transcripts,pval)
data <- head(results_transcripts, n=10)
write.csv(data, "../S288C_reference_genome_R64-2-1_20150113/significant.csv")
```

* Plot Distribution of transcript count per gene
```R
transcript_gene_table = indexes(bg_genome)$t2g
counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
temp <- hist(counts, breaks=50, col="bisque4", plot=FALSE)
plot(x = temp$mids, y = log10(temp$counts), type="h", xlab="Transcripts per gene", ylab="Log10 of Frequency", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)
```

* Plot Distribution of transcript lengths
```R
full_table <- texpr(bg_genome, 'all')
temp <- hist(full_table$length, breaks=50, col="steelblue", plot=FALSE)
plot(x = temp$mids, y = log10(temp$counts), type="h", xlab="Transcript length (bp)", ylab="Log10 of Frequency", main="Distribution of transcript lengths")
```

* Plot Distribution of transcript lengths
```R
full_table <- texpr(bg_genome, 'all')
temp <- hist(full_table$length, breaks=50, col="steelblue", plot=FALSE)
plot(x = temp$mids, y = log10(temp$counts), type="h", xlab="Transcript length (bp)", ylab="Log10 of Frequency", main="Distribution of transcript lengths")
```

* Get the gene expression and compute the maximum FPKM values for 4 libraries
```R
gene_expression = as.data.frame(gexpr(bg_genome_filt))
colnames(gene_expression) <- c("SRR1257637","SRR1257640","SRR1257793","SRR1259267")
max(gene_expression[,"SRR1257637"])
max(gene_expression[,"SRR1257640"])
max(gene_expression[,"SRR1257793"])
max(gene_expression[,"SRR1259267"])
```

* Plot Distribution of FPKM for all 4 libraries
```R
data_colors=c("tomato1","tomato2","wheat1","wheat2")

min_nonzero = 1
data_columns=c(1:4)
short_names=c("WT-1","WT-2","Isw2-1","Isw2-2")
boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 4 libraries")
```

* Plot expression values for pairs of wild-type replicates
```R
### First Pair
x = gene_expression[,"SRR1257637"]
y = gene_expression[,"SRR1257640"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (WT, Replicate 1)", ylab="FPKM (WT, Replicate 2)", main="Comparison of expression values for a pair of replicates")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

### Second Pair
x = gene_expression[,"SRR1257793"]
y = gene_expression[,"SRR1259267"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (Isw2, Replicate 1)", ylab="FPKM (Isw2, Replicate 2)", main="Comparison of expression values for a pair of replicates")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")
```

* Compute the Pearson correlation between all pairs of libraries; Plot Multidimensional Scaling (MDS) distance plots for the libraries

```R
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)
i = which(gene_expression[,"sum"] > 5)
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)
```

* View the distribution of differential expression values with significance p < 0.05
```R
bg_table = texpr(bg_genome_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) WT vs Isw2", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)
```

* Make plots of individual transcripts across samples
```R
ballgown::transcriptNames(bg_genome)[582] ## "MSTRG.278.4" 
plot(fpkm[12,] ~ pheno_data$population, border=c(1,2), main=paste(ballgown::transcriptNames(bg_genome)[582]),pch=19, xlab="Population", ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$population)), col=as.numeric(pheno_data$population))
```

* Display the grand expression values from UHR and HBR and mark those that are significantly differentially expressed
```R
gene_expression[,"WT"]=apply(gene_expression[,c(1, 2)], 1, mean)
gene_expression[,"Isw2"]=apply(gene_expression[,c(3, 4)], 1, mean)
x=log2(gene_expression[,"WT"]+min_nonzero)
y=log2(gene_expression[,"Isw2"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="WT FPKM (log2)", ylab="Isw2 FPKM (log2)", main="WT vs Isw2 FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

**Extract the annotation of each significant transcript using grep -a; Extract sequences of the most significant transcripts using Bedtools**  
```sh
grep -a ${SIGNIFICANT} ${REFERENCE}stringtie_merged.gtf > ${REFERENCE}significant1/${SIGNIFICANT}.gtf

bedtools getfasta -fi ${REFERENCE}S288C_reference_sequence_R64-2-1_20150113_new.fsa -bed ${REFERENCE}significant1/${SIGNIFICANT}.gtf -fo ${REFERENCE}significant1/${SIGNIFICANT}.fasta
```

## Code

```./scripts/bash_step1.sh```: Step 1 & 2 & 3  
```./scripts/bash_step2.sh```: Step 4 & 5  
```./scripts/bash_step3.sh```: Step 6  
```./scripts/R_script.sh```: Step 7  
(```./scripts/rscript1.R``` & ```./scripts/rscript2.R``` includes all R codes for Step 7)

***

# Conclusion
[Report](https://docs.google.com/document/d/1OVK1lC2Tv07apcZXxRsIEHGQw2ZwCAVIk3lZTvoO_bk/edit?usp=sharing)

**Detailed Directory Structure Listed Here:**  
```
.
|-- ballgown
|   |-- SRR1257637
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   |-- SRR1257637_S288C.gtf
|   |   `-- t_data.ctab
|   |-- SRR1257640
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   |-- SRR1257640_S288C.gtf
|   |   `-- t_data.ctab
|   |-- SRR1257793
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   |-- SRR1257793_S288C.gtf
|   |   `-- t_data.ctab
|   |-- SRR1259267
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   |-- SRR1259267_S288C.gtf
|   |   `-- t_data.ctab
|   |-- SRR1514794
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   |-- SRR1514794_S288C.gtf
|   |   `-- t_data.ctab
|   |-- SRR1514795
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   |-- SRR1514795_S288C.gtf
|   |   `-- t_data.ctab
|   |-- SRR1515155
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   |-- SRR1515155_S288C.gtf
|   |   `-- t_data.ctab
|   |-- SRR1515156
|   |   |-- e2t.ctab
|   |   |-- e_data.ctab
|   |   |-- i2t.ctab
|   |   |-- i_data.ctab
|   |   |-- SRR1515156_S288C.gtf
|   |   `-- t_data.ctab
|   |-- temp30
|   |   |-- SRR1257637
|   |   |   |-- e2t.ctab
|   |   |   |-- e_data.ctab
|   |   |   |-- i2t.ctab
|   |   |   |-- i_data.ctab
|   |   |   |-- SRR1257637_S288C.gtf
|   |   |   `-- t_data.ctab
|   |   |-- SRR1257640
|   |   |   |-- e2t.ctab
|   |   |   |-- e_data.ctab
|   |   |   |-- i2t.ctab
|   |   |   |-- i_data.ctab
|   |   |   |-- SRR1257640_S288C.gtf
|   |   |   `-- t_data.ctab
|   |   |-- SRR1257793
|   |   |   |-- e2t.ctab
|   |   |   |-- e_data.ctab
|   |   |   |-- i2t.ctab
|   |   |   |-- i_data.ctab
|   |   |   |-- SRR1257793_S288C.gtf
|   |   |   `-- t_data.ctab
|   |   `-- SRR1259267
|   |       |-- e2t.ctab
|   |       |-- e_data.ctab
|   |       |-- i2t.ctab
|   |       |-- i_data.ctab
|   |       |-- SRR1259267_S288C.gtf
|   |       `-- t_data.ctab
|   `-- temp37
|       |-- SRR1514794
|       |   |-- e2t.ctab
|       |   |-- e_data.ctab
|       |   |-- i2t.ctab
|       |   |-- i_data.ctab
|       |   |-- SRR1514794_S288C.gtf
|       |   `-- t_data.ctab
|       |-- SRR1514795
|       |   |-- e2t.ctab
|       |   |-- e_data.ctab
|       |   |-- i2t.ctab
|       |   |-- i_data.ctab
|       |   |-- SRR1514795_S288C.gtf
|       |   `-- t_data.ctab
|       |-- SRR1515155
|       |   |-- e2t.ctab
|       |   |-- e_data.ctab
|       |   |-- i2t.ctab
|       |   |-- i_data.ctab
|       |   |-- SRR1515155_S288C.gtf
|       |   `-- t_data.ctab
|       `-- SRR1515156
|           |-- e2t.ctab
|           |-- e_data.ctab
|           |-- i2t.ctab
|           |-- i_data.ctab
|           |-- SRR1515156_S288C.gtf
|           `-- t_data.ctab
|-- data
|   |-- samples.txt
|   |-- SRR1257637_1.fastq.gz
|   |-- SRR1257637_2.fastq.gz
|   |-- SRR1257640_1.fastq.gz
|   |-- SRR1257640_2.fastq.gz
|   |-- SRR1257793_1.fastq.gz
|   |-- SRR1257793_2.fastq.gz
|   |-- SRR1259267_1.fastq.gz
|   |-- SRR1259267_2.fastq.gz
|   |-- SRR1514794_1.fastq.gz
|   |-- SRR1514794_2.fastq.gz
|   |-- SRR1514795_1.fastq.gz
|   |-- SRR1514795_2.fastq.gz
|   |-- SRR1515155_1.fastq.gz
|   |-- SRR1515155_2.fastq.gz
|   |-- SRR1515156_1.fastq.gz
|   |-- SRR1515156_2.fastq.gz
|   |-- tri_copy
|   |   |-- SRR1257637_1_paired.fastq.gz
|   |   |-- SRR1257637_1_unpaired.fastq.gz
|   |   |-- SRR1257637_2_paired.fastq.gz
|   |   |-- SRR1257637_2_unpaired.fastq.gz
|   |   |-- SRR1257640_1_paired.fastq.gz
|   |   |-- SRR1257640_1_unpaired.fastq.gz
|   |   |-- SRR1257640_2_paired.fastq.gz
|   |   |-- SRR1257640_2_unpaired.fastq.gz
|   |   |-- SRR1257793_1_paired.fastq.gz
|   |   |-- SRR1257793_1_unpaired.fastq.gz
|   |   |-- SRR1257793_2_paired.fastq.gz
|   |   |-- SRR1257793_2_unpaired.fastq.gz
|   |   |-- SRR1259267_1_paired.fastq.gz
|   |   |-- SRR1259267_1_unpaired.fastq.gz
|   |   |-- SRR1259267_2_paired.fastq.gz
|   |   |-- SRR1259267_2_unpaired.fastq.gz
|   |   |-- SRR1514794_1_paired.fastq.gz
|   |   |-- SRR1514794_1_unpaired.fastq.gz
|   |   |-- SRR1514794_2_paired.fastq.gz
|   |   |-- SRR1514794_2_unpaired.fastq.gz
|   |   |-- SRR1514795_1_paired.fastq.gz
|   |   |-- SRR1514795_1_unpaired.fastq.gz
|   |   |-- SRR1514795_2_paired.fastq.gz
|   |   |-- SRR1514795_2_unpaired.fastq.gz
|   |   |-- SRR1515155_1_paired.fastq.gz
|   |   |-- SRR1515155_1_unpaired.fastq.gz
|   |   |-- SRR1515155_2_paired.fastq.gz
|   |   |-- SRR1515155_2_unpaired.fastq.gz
|   |   |-- SRR1515156_1_paired.fastq.gz
|   |   |-- SRR1515156_1_unpaired.fastq.gz
|   |   |-- SRR1515156_2_paired.fastq.gz
|   |   `-- SRR1515156_2_unpaired.fastq.gz
|   |-- trimmomatic
|   |   |-- SRR1257637_1_paired.fastq.gz
|   |   |-- SRR1257637_1_unpaired.fastq.gz
|   |   |-- SRR1257637_2_paired.fastq.gz
|   |   |-- SRR1257637_2_unpaired.fastq.gz
|   |   |-- SRR1257640_1_paired.fastq.gz
|   |   |-- SRR1257640_1_unpaired.fastq.gz
|   |   |-- SRR1257640_2_paired.fastq.gz
|   |   |-- SRR1257640_2_unpaired.fastq.gz
|   |   |-- SRR1257793_1_paired.fastq.gz
|   |   |-- SRR1257793_1_unpaired.fastq.gz
|   |   |-- SRR1257793_2_paired.fastq.gz
|   |   |-- SRR1257793_2_unpaired.fastq.gz
|   |   |-- SRR1259267_1_paired.fastq.gz
|   |   |-- SRR1259267_1_unpaired.fastq.gz
|   |   |-- SRR1259267_2_paired.fastq.gz
|   |   |-- SRR1259267_2_unpaired.fastq.gz
|   |   |-- SRR1514794_1_paired.fastq.gz
|   |   |-- SRR1514794_1_unpaired.fastq.gz
|   |   |-- SRR1514794_2_paired.fastq.gz
|   |   |-- SRR1514794_2_unpaired.fastq.gz
|   |   |-- SRR1514795_1_paired.fastq.gz
|   |   |-- SRR1514795_1_unpaired.fastq.gz
|   |   |-- SRR1514795_2_paired.fastq.gz
|   |   |-- SRR1514795_2_unpaired.fastq.gz
|   |   |-- SRR1515155_1_paired.fastq.gz
|   |   |-- SRR1515155_1_unpaired.fastq.gz
|   |   |-- SRR1515155_2_paired.fastq.gz
|   |   |-- SRR1515155_2_unpaired.fastq.gz
|   |   |-- SRR1515156_1_paired.fastq.gz
|   |   |-- SRR1515156_1_unpaired.fastq.gz
|   |   |-- SRR1515156_2_paired.fastq.gz
|   |   `-- SRR1515156_2_unpaired.fastq.gz
|   `-- Trimmomatic-0.36
|       |-- adapters
|       |   |-- NexteraPE-PE.fa
|       |   |-- TruSeq2-PE.fa
|       |   |-- TruSeq2-SE.fa
|       |   |-- TruSeq3-PE-2.fa
|       |   |-- TruSeq3-PE.fa
|       |   `-- TruSeq3-SE.fa
|       |-- LICENSE
|       `-- trimmomatic-0.36.jar
|-- results
|   |-- genome_gene_results1.csv
|   |-- genome_gene_results.csv
|   |-- genome_transcript_results1.csv
|   |-- genome_transcript_results.csv
|   |-- R_output.txt
|   |-- significant_transcripts.csv
|   |-- SRR1257637_S288C.bam
|   |-- SRR1257637_S288C.gtf
|   |-- SRR1257637_S288C.sam
|   |-- SRR1257640_S288C.bam
|   |-- SRR1257640_S288C.gtf
|   |-- SRR1257640_S288C.sam
|   |-- SRR1257793_S288C.bam
|   |-- SRR1257793_S288C.gtf
|   |-- SRR1257793_S288C.sam
|   |-- SRR1259267_S288C.bam
|   |-- SRR1259267_S288C.gtf
|   |-- SRR1259267_S288C.sam
|   |-- SRR1514794_S288C.bam
|   |-- SRR1514794_S288C.gtf
|   |-- SRR1514794_S288C.sam
|   |-- SRR1514795_S288C.bam
|   |-- SRR1514795_S288C.gtf
|   |-- SRR1514795_S288C.sam
|   |-- SRR1515155_S288C.bam
|   |-- SRR1515155_S288C.gtf
|   |-- SRR1515155_S288C.sam
|   |-- SRR1515156_S288C.bam
|   |-- SRR1515156_S288C.gtf
|   `-- SRR1515156_S288C.sam
|-- S288C_reference_genome_R64-2-1_20150113
|   |-- bedtools_output.fasta
|   |-- gene_association_R64-2-1_20150113.sgd
|   |-- grep.sh
|   |-- index
|   |   |-- S288C.1.ht2
|   |   |-- S288C.2.ht2
|   |   |-- S288C.3.ht2
|   |   |-- S288C.4.ht2
|   |   |-- S288C.5.ht2
|   |   |-- S288C.6.ht2
|   |   |-- S288C.7.ht2
|   |   `-- S288C.8.ht2
|   |-- merged.annotated.gtf
|   |-- merged.loci
|   |-- merged.stats
|   |-- merged.stringtie_merged.gtf.refmap
|   |-- merged.stringtie_merged.gtf.tmap
|   |-- merged.tracking
|   |-- mergelist.txt
|   |-- NotFeature_R64-2-1_20150113.fasta
|   |-- phenotype1.csv
|   |-- phenotype2.csv
|   |-- phenotype.csv
|   |-- rna_coding_R64-2-1_20150113.fasta
|   |-- S288C.gtf
|   |-- S288C_reference_sequence_R64-2-1_20150113
|   |-- S288C_reference_sequence_R64-2-1_20150113.fsa
|   |-- S288C_reference_sequence_R64-2-1_20150113.fsa.fai
|   |-- S288C_reference_sequence_R64-2-1_20150113_new.fsa
|   |-- S288C_reference_sequence_R64-2-1_20150113_new.fsa.fai
|   |-- saccharomyces_cerevisiae_R64-2-1_20150113.gff
|   |-- sed
|   |-- significant1
|   |   |-- MSTRG.107.1.fasta
|   |   |-- MSTRG.107.1.gtf
|   |   |-- MSTRG.179.1.fasta
|   |   |-- MSTRG.179.1.gtf
|   |   |-- MSTRG.196.1.fasta
|   |   |-- MSTRG.196.1.gtf
|   |   |-- MSTRG.258.1.fasta
|   |   |-- MSTRG.258.1.gtf
|   |   |-- MSTRG.278.4.fasta
|   |   |-- MSTRG.278.4.gtf
|   |   |-- MSTRG.30.2.fasta
|   |   |-- MSTRG.30.2.gtf
|   |   |-- MSTRG.344.1.fasta
|   |   |-- MSTRG.344.1.gtf
|   |   |-- MSTRG.66.13.fasta
|   |   |-- MSTRG.66.13.gtf
|   |   |-- MSTRG.8.3.fasta
|   |   |-- MSTRG.8.3.gtf
|   |   |-- MSTRG.88.1.fasta
|   |   `-- MSTRG.88.1.gtf
|   |-- significant1.txt
|   |-- significant2
|   |   |-- MSTRG.129.1.fasta
|   |   |-- MSTRG.129.1.gtf
|   |   |-- MSTRG.183.3.fasta
|   |   |-- MSTRG.183.3.gtf
|   |   |-- MSTRG.197.2.fasta
|   |   |-- MSTRG.197.2.gtf
|   |   |-- MSTRG.304.15.fasta
|   |   |-- MSTRG.304.15.gtf
|   |   |-- MSTRG.315.1.fasta
|   |   |-- MSTRG.315.1.gtf
|   |   |-- MSTRG.376.1.fasta
|   |   |-- MSTRG.376.1.gtf
|   |   |-- MSTRG.395.1.fasta
|   |   |-- MSTRG.395.1.gtf
|   |   |-- MSTRG.42.5.fasta
|   |   |-- MSTRG.42.5.gtf
|   |   |-- MSTRG.66.26.fasta
|   |   |-- MSTRG.66.26.gtf
|   |   |-- MSTRG.92.1.fasta
|   |   `-- MSTRG.92.1.gtf
|   |-- significant2.txt
|   |-- stringtie_merged.gtf
|   `-- tmp.XXFy1XkF
|-- scripts
|   |-- bash_step1.sh
|   |-- bash_step2.sh
|   |-- bash_step3.sh
|   |-- bash_step4-1.sh
|   |-- bash_step4-2.sh
|   |-- clean_fsa.sh
|   |-- download_fastq.sh
|   |-- fastqc.sh
|   |-- gffread.sh
|   |-- index_builder.sh
|   |-- rscript1.R
|   |-- rscript2.R
|   |-- R_script.sh
|   `-- trimmomatic.sh
`-- Trimmomatic-0.36
    |-- adapters
    |   |-- NexteraPE-PE.fa
    |   |-- TruSeq2-PE.fa
    |   |-- TruSeq2-SE.fa
    |   |-- TruSeq3-PE-2.fa
    |   |-- TruSeq3-PE.fa
    |   `-- TruSeq3-SE.fa
    |-- LICENSE
    `-- trimmomatic-0.36.jar
```
