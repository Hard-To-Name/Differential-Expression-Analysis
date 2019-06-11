library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

###8. Load phenotype data for the samples
pheno_data = read.csv("../S288C_reference_genome_R64-2-1_20150113/phenotype2.csv")

###9. Read in the expression data that were calculated by StringTie
bg_genome = ballgown(dataDir = "../ballgown/temp37", samplePattern = "SRR", pData = pheno_data) 

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
head(results_transcripts, n=10) ## select top 10 transcripts with the smallest p values

### Plot Distribution of transcript count per gene
transcript_gene_table = indexes(bg_genome)$t2g
counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
temp <- hist(counts, breaks=50, col="bisque4", plot=FALSE)
plot(x = temp$mids, y = log10(temp$counts), type="h", xlab="Transcripts per gene", ylab="Log10 of Frequency", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)

### Plot Distribution of transcript lengths
full_table <- texpr(bg_genome, 'all')
temp <- hist(full_table$length, breaks=50, col="steelblue", plot=FALSE)
plot(x = temp$mids, y = log10(temp$counts), type="h", xlab="Transcript length (bp)", ylab="Log10 of Frequency", main="Distribution of transcript lengths")

### FPKM values
gene_expression = as.data.frame(gexpr(bg_genome_filt))
colnames(gene_expression) <- c("SRR1514794","SRR1514795","SRR1515155","SRR1515156")
max(gene_expression[,"SRR1514794"])
max(gene_expression[,"SRR1514795"])
max(gene_expression[,"SRR1515155"])
max(gene_expression[,"SRR1515156"])

### Plot Distribution of FPKM for all 4 libraries
data_colors=c("tomato1","tomato2","wheat1","wheat2")

min_nonzero = 1
data_columns=c(1:4)
short_names=c("Rsc-1","Rsc-2","Ino80-1","Ino80-2")
boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 4 libraries")

### Plot expression values for a pair of Rsc mutant replicates
x = gene_expression[,"SRR1514794"]
y = gene_expression[,"SRR1514795"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (Rsc, Replicate 1)", ylab="FPKM (Rsc, Replicate 2)", main="Comparison of expression values for a pair of replicates")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

### Plot expression values for a pair of Ino80 mutant replicates
x = gene_expression[,"SRR1515155"]
y = gene_expression[,"SRR1515156"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (Ino80, Replicate 1)", ylab="FPKM (Ino80, Replicate 2)", main="Comparison of expression values for a pair of replicates")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")


### MDS distance plot
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)
i = which(gene_expression[,"sum"] > 5)
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

### View the distribution of differential expression values with significance p < 0.1
bg_table = texpr(bg_genome_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

sig=which(results_genes$pval<0.1)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) Rsc vs Ino80", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

###19. Make plots of individual transcripts across samples
ballgown::transcriptNames(bg_genome)[582] ## "MSTRG.278.4" 
plot(fpkm[12,] ~ pheno_data$population, border=c(1,2), main=paste(ballgown::transcriptNames(bg_genome)[582]),pch=19, xlab="Population", ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$population)), col=as.numeric(pheno_data$population))

### Display the grand expression values from UHR and HBR and mark those that are significantly differentially expressed
gene_expression[,"Rsc"]=apply(gene_expression[,c(1, 2)], 1, mean)
gene_expression[,"Ino80"]=apply(gene_expression[,c(3, 4)], 1, mean)
x=log2(gene_expression[,"Rsc"]+min_nonzero)
y=log2(gene_expression[,"Ino80"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="WT FPKM (log2)", ylab="Lsw2 FPKM (log2)", main="Rsc vs Ino80 FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

### significant transcripts with a log2 fold-change >= 2
sigpi = which(results_genes[,"pval"]<0.1)
sigp = results_genes[sigpi,]
sigde = which(abs(sigp[,"de"]) >= 1)
sig_tn_de = sigp[sigde,]

o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)
output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
#write.table(output, file="SigDE.txt", sep="\t", row.names=FALSE, quote=FALSE)
#View selected columns of the first 25 lines of output
output[1:25,c(1,4,5)]
