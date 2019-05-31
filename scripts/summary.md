## Methods

**Data**

The release 64-2-1 of Saccharomyces cerevisiae S228C reference genome is obtained from Saccharomyces genome database ([https://downloads.yeastgenome.org/sequence/S288C\_reference/genome\_releases/](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/)). The strain is selected to be non-flocculent with a minimal set of nutritional requirements. The Illumina data of all 8 samples is obtained from European Bioinformatic Institute (EBI) and can be downloaded by script available at ([https://github.com/Hard-To-Name/Differential-Expression-Analysis/blob/master/scripts/download\_fastq.sh](https://github.com/Hard-To-Name/Differential-Expression-Analysis/blob/master/scripts/download_fastq.sh)).

**Quality Measurement of Sample Sequence**

Statistics of all 8 samples were measured by FASTQC (Andrews, S. 2014) version 0.11.8 including adapter content, duplication levels, per base n content, per base sequence quality, per base sequence content, per sequence gc content, per sequence quality, per tile quality, and sequence length distribution. The quality statistics were used to determine characteristics of data and to determine the parameters for data trimming.

**Data Trimming**

Trimmomatic (Bolger et al. 2014) version 0.39 was used to trim the 8x illumina paired-end data of all 8 samples from Sequence Read Archive (SRA) by removing adapters, removing leading and trailing bases with quality below 0, scanning the read with a 4-base wide sliding window with cutting when the average quality per base drops below 15, dropping reads below the 20 bases long. HISAT2 (Kim et al. 2015) version 2.1.0 was used to build HISAT2 indices for the reference genome.

**Genome Assemblies**

Hisat2 was used to align the 8x Illumina RNA-seq paired-end reads to the reference genome, resulting in sequence alignment map (.sam) files. Samtools (Li et al. 2009) version 1.9 was used to convert all .sam files to binary alignment map (.bam) files. StringTie (Pertea et al. 2015) version 1.3.6 was used to assemble RNA-Seq alignments into potential transcripts and output an estimated gene transfer format (.gtf) file for each sample, merge all transcripts generated into a non-redundant set of transcripts and output a merged .gtf file. Gffcompare under Cufflinks ([http://ccb.jhu.edu/software/stringtie/gffcompare.shtml](http://ccb.jhu.edu/software/stringtie/gffcompare.shtml)) version 0.11.2 was used to compare and compute the estimated accuracy of the merged .gtf file and .gtf files of each sample, to check how the predicted transcripts relate to an annotation file. Finally StringTie was used to estimate transcript abundances, create table counts, generate .gtf file and .bam file for each sample based on the merged annotation file.

**Gene Expression Analysis**

All 8 samples were separated into two groups: one group has four samples under temperature 30 centigrades, with two samples belonging to population Lsw2 and the other two samples belonging to population WT; the other group has four samples under temperature 37 centigrades, with two samples belonging to population Rsc and the other two samples belonging to population Ino80.

Phenotype data and the expression data calculated by StringTie were read into R. Several steps were performed using R package Ballgown (Frazee, A.C. et al. 2015) version 2.10.0 to prepare the data: first, filter out low-abundance genes which have a variance among all samples less than 1; second, identify transcripts and genes that show statistically significant differences with population as the covariate of interest under the same temperature and acquire the fold changes of the identified transcripts and genes with FPKM as measurement; third, arrange both the identified transcripts and the identified genes by increasing p-values; fourth, save the resulting transcripts and genes to separate .csv files.

Then transcripts and genes with q-values less than 0.05 were identified and considered as having significant statistical difference. Ballgown was used to examine the distribution of gene abundances in FPKM across all samples. The distribution of expression data was then normalized by log transformation and plotted.

After that transcripts were plotted individually. The structure and expression levels of all transcripts that share the same gene locus within sample SRR1257637 were plotted. And the average expression levels for all transcript of gene MSTRG.128 was plotted.



## Results

**Data Quality**

**Data Trimming Results**

**Genome Assembly**

**Differently Expressed Transcripts**

**Differently Expressed Genes**

**Gene Expression Analysis**



## Conclusion



## Work Cited

Bolger, A. M., Lohse, M., &amp; Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. [https://doi.org/10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170)

Engel, S. R., Dietrich, F. S., Fisk, D. G., Binkley, G., Balakrishnan, R., Costanzo, M. C., … Cherry, J. M. (2014). The Reference Genome Sequence of Saccharomyces cerevisiae : Then and Now. G3&amp;amp;#58; Genes|Genomes|Genetics, 4(3), 389–398. [https://doi.org/10.1534/g3.113.008995](https://doi.org/10.1534/g3.113.008995)

Frazee, A. C., Pertea, G., Jaffe, A. E., Langmead, B., Salzberg, S. L., &amp; Leek, J. T. (2015). Ballgown bridges the gap between transcriptome assembly and expression analysis. Nature Biotechnology, 33(3), 243–246. [https://doi.org/10.1038/nbt.3172](https://doi.org/10.1038/nbt.3172)

Kim, D., Langmead, B., &amp; Salzberg, S. L. (2015). HISAT: A fast spliced aligner with low memory requirements. Nature Methods, 12(4), 357–360. [https://doi.org/10.1038/nmeth.3317](https://doi.org/10.1038/nmeth.3317)

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. [https://doi.org/10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)

Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., &amp; Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. [https://doi.org/10.1038/nprot.2016.095](https://doi.org/10.1038/nprot.2016.095)

Pertea, M., Pertea, G. M., Antonescu, C. M., Chang, T. C., Mendell, J. T., &amp; Salzberg, S. L. (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nature Biotechnology, 33(3), 290–295. [https://doi.org/10.1038/nbt.3122](https://doi.org/10.1038/nbt.3122)

Trapnell, C., Roberts, A., Goff, L., Pertea, G., Kim, D., Kelley, D. R., … Pachter, L. (2014, January 1). Erratum: Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks (Nature Protocols (2012) 7 (562-578)). Nature Protocols, Vol. 9, p. 2513. [https://doi.org/10.1038/nprot.2012.016](https://doi.org/10.1038/nprot.2012.016)

Wingett, S. W., &amp; Andrews, S. (2018). FastQ Screen: A tool for multi-genome mapping and quality control. F1000Research, 7, 1338. [https://doi.org/10.12688/f1000research.15931.2](https://doi.org/10.12688/f1000research.15931.2)

Andrews, S. (2014). FastQC A Quality Control tool for High Throughput Sequence Data. [http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.)
