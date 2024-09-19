# ITS analysis of 'PLUS-ONE' project
# DADA2 pipeline
# Sequencing by SeqCoast, Portsmouth, NH
# V3-V4 on a NextSeq

# Seagrass Microbiome Analyses - ITS Only
# "Plus-One" Project
# Sequencing by Seqcoast 
# Updated: August 2024 by ANB

# Prepare your environment----
# Libraries
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ShortRead) 

# Check package versions
packageVersion("dada2")
packageVersion("ShortRead")
packageVersion("Biostrings")

# To install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead")

# Switch when processing ITS
setwd("~/Dropbox/SEAGRASS-MICROBIOME/SEAGRASS-SEQUENCES/RAW-DATA/ITS-DATA/")
path <- "~/Dropbox/SEAGRASS-MICROBIOME/SEAGRASS-SEQUENCES/RAW-DATA/ITS-DATA"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

# Identify primers
FWD <- "GCATCGATGAAGAACGCAGC"  ## CHANGE ME to your forward primer sequence
REV <- "TCCTCCGCTTATTGATATGC"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients,primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits,fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Remove the primers using cutadapt =============================================
# Use conda env list to figure out the paths
# Make sure it's pointing to the actual program and not just the directory
cutadapt <- "/Users/ashleybulseco/miniconda3/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R, it should come back with version
# If not, then it isn't able to find where your cutadapt version is

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
# TO remove 0-length sequences (https://github.com/benjjneb/dada2/issues/1595)
# for(i in seq_along(fnFs)) {
#   system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
#                              "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
#                              fnFs.filtN[i], fnRs.filtN[i])) # input files
# }

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-m", 20, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Check that it worked
# This should be 0's across the board, indicating that no primers were detected
# Success!
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients,primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits,fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Then you can proceed with the rest of the DADA2 pipeline
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Inspect quality profiles
plotQualityProfile(cutFs[1:2])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)
saveRDS(out,"ITS-out.rds")
out <- readRDS("ITS-out.rds")

# Learn the error rates
ITS_errF <- learnErrors(filtFs, multithread=TRUE)
ITS_errR <- learnErrors(filtRs, multithread = TRUE)
saveRDS(ITS_errF,"ITS_errF.rds")
saveRDS(ITS_errR,"ITS_errR.rds")
plotErrors(ITS_errF, nominalQ=TRUE)
plotErrors(ITS_errR, nominalQ=TRUE)
ITS_errF <- readRDS("ITS_errF.rds")
ITS_errR <- readRDS("ITS_errR.rds")

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
saveRDS(derepFs, "derepFS.rds")
saveRDS(derepRs, "derepRs.rds")
derepFs <- readRDS("derepFS.rds")
derepRs <- readRDS("derepRS.rds")

# We are now ready to apply the core sample inference algorithm to the dereplicated data.
ITS_dadaFs <- dada(derepFs, err=ITS_errF, multithread=TRUE)
ITS_dadaRs <- dada(derepRs, err=ITS_errR, multithread=TRUE)
saveRDS(ITS_dadaFs, "ITS_dadaFs")
saveRDS(ITS_dadaRs, "ITS_dadaRs")
ITS_dadaFs <- readRDS("ITS_dadaFs")
ITS_dadaRs <- readRDS("ITS_dadaRs")

# Merge reads
ITS_mergers <- mergePairs(ITS_dadaFs, filtFs, ITS_dadaRs, filtRs, verbose=TRUE)
head(ITS_mergers[[1]])
saveRDS(ITS_mergers, "ITS_mergers.rds")
# To read back in
ITS_mergers <- readRDS("ITS_mergers.rds")

# Sequence table
ITS_seqtab <- makeSequenceTable(ITS_mergers)
dim(ITS_seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(ITS_seqtab)))

# Detect and remove chimeras
ITS_seqtab.nochim <- removeBimeraDenovo(ITS_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(ITS_seqtab.nochim)))
saveRDS(ITS_seqtab.nochim, "ITS_seqtab.nochim.rds")
ITS_seqtab.nochim <- readRDS("ITS_seqtab.nochim.rds")

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(ITS_dadaFs, getN), sapply(ITS_dadaRs, getN), sapply(ITS_mergers, getN),
               rowSums(ITS_seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign Taxonomy
# https://benjjneb.github.io/dada2/training.html
# https://unite.ut.ee/repository.php
# Download from here: https://doi.plutof.ut.ee/doi/10.15156/BIO/2959332
unite.ref <- "~/Dropbox/SEAGRASS-MICROBIOME/SEAGRASS-SEQUENCES/RAW-DATA/ITS-DATA/sh_general_release_dynamic_all_04.04.2024.fasta"  # CHANGE ME to location on your machine
ITS_taxa <- assignTaxonomy(ITS_seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
saveRDS(ITS_taxa, "ITS_taxa.rds")
ITS_taxa <- readRDS("ITS_taxa.rds")

# Inspect
ITS_taxa.print <- ITS_taxa  # Removing sequence rownames for display only
rownames(ITS_taxa.print) <- NULL
ITS_taxa.print

# Export to better understand what's going on
# tax table:
ITS_asv_tax <- ITS_taxa
row.names(ITS_asv_tax) <- sub(">", "", asv_headers)
write.table(ITS_asv_tax, "ITS_ASVs_taxonomy_091624.tsv", sep="\t", quote=F, col.names=NA)

# # Load necessary libraries
# library(dada2)
# 
# # Assume your taxonomy table is called 'taxa'
# # taxa <- assignTaxonomy(seqtab, "path/to/reference_database.fa.gz")
# 
# # If your taxonomy table is a matrix, convert it to a data frame
# taxa_df <- as.data.frame(ITS_taxa)
# 
# # Export the taxonomy table to a text file
# write.table(taxa_df, file = "ITS_taxa_assignment.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
# 
# # Getting really poor taxonomic assignment, very few sequences assigned below Fungi

# Exporting ITS for Phyloseq
## Preparing for Phyloseq----------------
# Constructing a sample table that we can use for Phyloseq import
seqs <- getSequences(ITS_seqtab.nochim)
asv_seqs <- colnames(ITS_seqtab.nochim)
asv_headers <- vector(dim(ITS_seqtab.nochim)[2], mode="character")

for (i in 1:dim(ITS_seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(ITS_seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- ITS_taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

# Sample names:
write.table(sample.names, "sample_names.txt", sep="\t", quote=F, col.names=NA)
