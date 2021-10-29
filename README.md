# Hormone-Microbiome-DADA2-Data-Analysis
Human Steroid Hormones Influence Gut Microbial Community Structure - DADA2 Analysis

**Software Required**:
R/RStudio

DADA2: for Illumina-sequenced paired-end fastq files where barcodes/adapters have been removed.
[DADA2 Pipeline](https://benjjneb.github.io/dada2/tutorial_1_8.html)

# To install DADA2
**Need to first install Bioconductor and then dada2**:

`if (!requireNamespace("BiocManager", quietly = TRUE))  
    install.packages("BiocManager")`  
`BiocManager::install("dada2")`

Documentation: https://www.bioconductor.org/packages/release/bioc/html/dada2.html

Additional packages required for downstream analysis: phyloseq and ggplot2
