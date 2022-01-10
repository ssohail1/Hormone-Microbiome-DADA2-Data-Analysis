# Hormone-Microbiome-DADA2-Data-Analysis
Human Steroid Hormones Influence Gut Microbial Community Structure - DADA2 Analysis

**Software/Packages Required**:
- R: Bioconductor, dada2, ggplot2, phyloseq

### Installing DADA2
DADA2: for Illumina-sequenced paired-end fastq files where barcodes/adapters have been removed.  
[DADA2 Pipeline](https://benjjneb.github.io/dada2/tutorial_1_8.html)

**Need to first install Bioconductor and then dada2**:

`if (!requireNamespace("BiocManager", quietly = TRUE))`  
`install.packages("BiocManager")`  
`BiocManager::install("dada2")`

Documentation: https://www.bioconductor.org/packages/release/bioc/html/dada2.html

### Files in repo
data folder:
  - HormoneRscript_upd.R: R script with the [DADA2](https://benjjneb.github.io/dada2/tutorial_1_8.html) and phylogenetic analyses
  - filtfiles folder: contains the data used in the DADA2 and phylogenetic analyses.

### How to use
Clone repository into personal directory using this command,  
`git clone https://github.com/ssohail1/Hormone-Microbiome-DADA2-Data-Analysis.git`

To move into bc_microbiome-project directory use `cd`,  
`cd Hormone-Microbiome-DADA2-Data-Analysis`
