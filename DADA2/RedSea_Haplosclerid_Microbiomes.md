# R-Script DADA2 Pipeline


### Setup
```python
knitr::opts_chunk$set(echo = TRUE)
```

### Preparation
```python
if(!requireNamespace("dada2", quietly = TRUE)){
  print("dada2 not installed, installing now...")
  if(!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    install_exit_status <- BiocManager::install(c("dada2", "phyloseq"), lib="/PATH/TO/FOLDER/R/x86_64-suse-linux-gnu-library/4.3")
  }else{
    install_exit_status <- BiocManager::install(c("dada2","phyloseq"), lib="/PATH/TO/FOLDER/R/x86_64-suse-linux-gnu-library/4.3")  
  }
}else{
  print("dada2 there, loading...")
  library(dada2)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(egg)
  library(ape)
  library(RColorBrewer)
  library(phytools)
}
```
### Load the clean reads

Samples "GW3237", "GW3325", "GW3404", "GW3416", "GW3494", "GW3580", "GW4166", "GW4251", "GW5875", "GW6059", "GW6104", "GW6140", "GW6158", and "GW6196" were discarded due to low read counts (<5000).

Samples "GW3524", "GW4232", "GW5899", and "GW6180" were discarded due to missing UCE data.

```python
path <- "./Data/Red_Sea_Prokaryotic_Microbiomes/Reads/"
#list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.clean.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.clean.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

### Filter reads

Place filtered files in filtered/ subdirectory
```python
filtFs <- file.path(path, "Dada2_Filtered", paste0(sample.names, "_F_dada2.clean.fastq.gz"))
filtRs <- file.path(path, "Dada2_Filtered", paste0(sample.names, "_R_dada2.clean.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

filter_output <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
```

Learn Errors
```python
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
```

Sample inference
```python
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

Merge reads
```python
merged_pairs <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

Get sequence table
```python
seqtab <- makeSequenceTable(merged_pairs)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)
```

Track reads
```python
getN <- function(x) sum(getUniques(x))
track <- cbind(filter_output, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged_pairs, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim)/filter_output[,2])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "surviving reads")
rownames(track) <- sample.names
#rownames(track[track[,1]<5000,])
```

### Assign Taxonomy
```python
red_sea_taxa <- assignTaxonomy(seqtab.nochim,"./Data/Red_Sea_Prokaryotic_Microbiomes/SilvaDB/v138/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

red_sea_taxa <- addSpecies(red_sea_taxa, "./Data/Red_Sea_Prokaryotic_Microbiomes/SilvaDB/v138/silva_species_assignment_v138.1.fa.gz")
```

Inspect taxonomy table
```python
red_sea_taxa.print <- red_sea_taxa # Removing sequence rownames for display only
rownames(red_sea_taxa.print) <- NULL
head(red_sea_taxa.print)
```

Extract ASV table
```python
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "./Data/ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "./Data/ASVs_Counts.tsv", sep="\t", quote=F, col.names=NA)

#taxonomy table:
row.names(red_sea_taxa.print) <- sub(">", "", asv_headers)
write.table(red_sea_taxa.print, "./Data/ASVs_Taxonomy.tsv", sep = "\t", quote=F, col.names=NA)


# #filter mitochondria, column 5 contains the family info
mt_ASVs <- sub(">", "", rownames(red_sea_taxa.print[which(red_sea_taxa.print[,5] == "Mitochondria"),]))
chloro_ASVs <- sub(">", "", rownames(red_sea_taxa.print[which(red_sea_taxa.print[,4] == "Chloroplast"),]))

asv_tab_noMt_noChloro <- asv_tab[!(row.names(asv_tab) %in% c(mt_ASVs, chloro_ASVs)),]

# count table:
write.table(asv_tab_noMt_noChloro, "./Data/ASVs_Counts_noMt_noChloro.tsv", sep="\t", quote=F, col.names=NA)

#taxonomy table:
red_sea_taxa.print_noMt_noChloro <- red_sea_taxa.print[!(row.names(red_sea_taxa.print) %in% c(mt_ASVs, chloro_ASVs)),]

write.table(red_sea_taxa.print_noMt_noChloro, "./Data/ASVs_Taxonomy_noMt_noChloro.tsv", sep = "\t", quote=F, col.names=NA)
```

