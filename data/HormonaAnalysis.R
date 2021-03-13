# Set the working directory
library(dada2)
setwd("/media/mbb/datapartition/Sidra/HormoneMicrobioTestingwithDada2/PrimerHormoneDADA") #change to directory of your files
dir()
path <- "/media/mbb/datapartition/Sidra/HormoneMicrobioTestingwithDada2/PrimerHormoneDADA"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001_primer.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001_primer.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# To get the quality profile of the forward reads:
plotQualityProfile(fnFs[1:2])

# To get the quality profile of the reverse reads:
plotQualityProfile(fnRs[1:2])

# Assigning filenames for the filtered fastq.gz files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,100),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

# Learn the Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

# Applying the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)


dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


dadaFs[[1]]
## dada-class: object describing DADA2 denoising results
## 216 sequence variants were inferred from 13523 input unique sequences.
## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Merging paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

# Constructing sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
## [1]   212 20897

table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE, 
                                    +                                     minFoldParentOverAbundance= 30)
## minFoldParentOverAbundance= 30 added so to avoid removing too many sequences from input data
## Identified 8964 bimeras out of 20897 input sequences.

sum(seqtab.nochim)/sum(seqtab)
## [1] 0.9593827 <-- after minFoldParentOverAbundance

# Track reads thru pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/media/mbb/datapartition/Sidra/silva_nr_v132_train_set.fa.gz", multithread = TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

### Skipped this part b/c no mock sample in hormone data
# Evaluate accuracy
#unqs.mock <- seqtab.nochim["Mock",]
#unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
#cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
#mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
#match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
#cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

# Phyloseq
## Need to download phyloseq first
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phyloseq") 

## Then load library
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
# this is for making the metadata,  I have my metadata available below
theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
Sample_Name <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(Sample_Name,1,1)
subject <- substr(Sample_Name,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(sample_name=Sample_Name, Name=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

# setting the metadata so to compare two variables
samdf <- read.table("~/Desktop/TestData.txt", header=TRUE)

# complete Metadata
samdf <- read.table("~/Desktop/Metadata.txt", header=TRUE)

#for comparing the Estradiol samples
samdf2 <- read.table("~/Desktop/Estradiolsamples.txt", header=TRUE)


#Phyloseq analysis for getting abundance bar plots
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
#Shannon-Simpson plot Alpha-diversity plot
plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"), color="Timepoint")
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray", color="Timepoint")
#Bray NMDS Beta-diversity plot
plot_ordination(ps.prop, ord.nmds.bray, color="Timepoint", title="Bray NMDS")
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
#Can have the fill = be at the Phylum or Genus level as well
plot_bar(ps.top20, x="Treatment", fill="Family") + facet_wrap(~Timepoint, scales="free_x")
# End of dada2 Tutorial

# RAREFACTION
## Do rarecurve command first then rarefy command 
install.packages("vegan")
library(vegan)
data(seqtab.nochim)
S <- specnumber(seqtab.nochim) # observed number of species
(raremax <- min(rowSums(seqtab.nochim)))
plot(sort(rowSums(seqtab.nochim))[1:34])
Srare <- rarefy(seqtab.nochim, raremax)
plot(S, Srare, xlab = "Observed Number of Species", ylab = "Rarefied No. of Species")
abline(0, 0.1117482900000299990)
#Get the rarecurve at this step
rarecurve(seqtab.nochim, step = 100, sample = raremax, col = "blue", cex = 0.6)
S
rarefied_data<-rrarefy(seqtab.nochim, 7000)
