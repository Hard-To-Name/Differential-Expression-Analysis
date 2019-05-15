# File Structure
  
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
└── S288C_reference_genome_R64-2-1_20150113  
|   ├── index  
|   |   └── *.ht2  
|   ├── .gtf of reference genome  
|   ├── .fsa of reference genome before and after cleaning headers  
|   ├── phenotype.csv for ballgown  
|   └── merge information  
├── results  
|   ├── *.sam before merge  
|   ├── *.bam before merge  
|   └── *.gtf before merge  
└── ballgown  
    ├── *.bam after merge  
    └── *.gtf after merge  
