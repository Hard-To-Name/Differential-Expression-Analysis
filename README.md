# Differential-Expression-Analysis

***

## Pipeline

**1. Map the reads for each sample to the reference genome**  
  
**2. Sort and convert the SAM files to BAM**  

**3. Assemble transcripts for each sample**  

**4. Merge transcripts from all samples**  

**5. Examine how the transcripts compare with the reference annotation (optional)**  

**6. Estimate transcript abundances and create table counts for Ballgown**  

**7. Load relevant R packages**  

**8. Load the phenotype data for the samples**  

**9. Read in te expression data that were calculated by StringTie.**  

**10. Filter to remove low-abundance genes**  

**11. Identify transcripts that show statistically significant differences between groups**  

**12. Identify genes that show statistically significant differences between groups**  

**13. Add gene names and gene IDs to the ```results_transcripts``` data frame**  

**14. Sort the results from the smallest P value to the largest**  

**15. Write the results to a ```csv``` file that can be shared and distributed**  

**16. Identify transcripts and genes with a q value < 0.05**  

**17. Plotting (optional)**  

**18. Show the distribution of gene abundances (measured as FPKM values) across samples**  

**19. Make plots of individual transcripts across samples**  

**20. Plot the structure and expression levels in a sample of all transcripts that share the same gene locus**  

**21. Plot the average expression levels for all transcripts of a gene within different groups using the ```plotMeans```function**
