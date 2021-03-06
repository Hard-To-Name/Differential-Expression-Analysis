#!/usr/bin/env Rscript

#$PHENOTYPEDATA="csv_of_phenotype_data_information.csv"
#$SAMPLENAMEPATTERN="common_sample_name_pattern"
#$FILTERINGLINE="ex_rowVars(texpr(bg_chrX)) >1"

###7. Load relevant R packages
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

###8. Load phenotype data for the samples
pheno_data = read.csv("../S288C_reference_genome_R64-2-1_20150113/phenotype1.csv")

###9. Read in the expression data that were calculated by StringTie
bg_genome = ballgown(dataDir = "../ballgown/temp30", samplePattern = "SRR", pData = pheno_data) 

###10. Filter to remove low-abundance genes
bg_genome_filt = subset(bg_genome, "rowVars(texpr(bg_genome)) > 1", genomesubset=TRUE)

###11. Identify transcripts that show statistically significant differences between groups
results_transcripts = stattest(bg_genome_filt, feature="transcript",covariate="population", getFC=TRUE, meas="FPKM")

###12. Identify genes that show statistically significant differences between groups
results_genes = stattest(bg_genome_filt, feature="gene", covariate="population", getFC=TRUE, meas="FPKM")

###13. Add gene names and gene IDs to the results_transcripts data frame
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_genome_filt), geneIDs=ballgown::geneIDs(bg_genome_filt), results_transcripts)

###14. Sort the results from the smallest P value to the largest
results_transcripts = arrange(results_transcripts,pval) 
results_genes = arrange(results_genes,pval)

###15. Write the results to a csv file that can be shared and distributed
write.csv(results_transcripts, "../results/genome_transcript_results1.csv", row.names=FALSE) 
write.csv(results_genes, "../results/genome_gene_results1.csv", row.names=FALSE) 

###16. Identify transcripts and genes with a q value < 0.05
subset(results_transcripts,results_transcripts$qval<0.05) 
subset(results_genes,results_genes$qval<0.05)

###17. Plotting (optional)
tropical= c('darkorange', 'dodgerblue', 'hotpink', 'yellow') 
palette(tropical)

###18. Show the distribution of gene abundances (measured as FPKM values) across samples
fpkm = texpr(bg_genome,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm,col=as.numeric(pheno_data$population),las=2,ylab='log2(FPKM+1)')

###19. Make plots of individual transcripts across samples
ballgown::transcriptNames(bg_genome)[582] ##      1200 ## "MSTRG.278.4" 
ballgown::geneNames(bg_genome)[582] ##      1200 ## "." 
plot(fpkm[12,] ~ pheno_data$population, border=c(1,2), main=paste(ballgown::transcriptNames(bg_genome)[582]),pch=19, xlab="Population", ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$population)), col=as.numeric(pheno_data$population))

###20. Plot the structure and expression levels in a sample of all transcripts that share the same gene locus
plotTranscripts(ballgown::geneIDs(bg_genome)[1729], bg_genome, main=c('Gene XIST in sample SRR1257637'), sample=c('SRR1257640'))

###21. Plot the average expression levels for all transcripts of a gene within different groups using the plotMeansfunction
plotMeans('MSTRG.157', bg_genome_filt,groupvar="population",legend=FALSE) ## qval = 0.04
