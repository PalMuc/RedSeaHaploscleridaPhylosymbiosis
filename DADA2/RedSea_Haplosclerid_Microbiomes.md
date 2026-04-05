# R-Script DADA2 Pipeline

Please check this tutorial before you proceed: https://benjjneb.github.io/dada2/tutorial.html
and: https://www.bioconductor.org/packages//release/bioc/vignettes/dada2/inst/doc/dada2-intro.html

### Setup
```r
knitr::opts_chunk$set(echo = TRUE)
```

### Preparation
```r
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
```r
path <- "./Data/Red_Sea_Prokaryotic_Microbiomes/Reads/"

fnFs <- sort(list.files(path, pattern="_R1_001.clean.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.clean.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

### Filter reads

Place filtered files in filtered/ subdirectory.
```r
filtFs <- file.path(path, "Dada2_Filtered", paste0(sample.names, "_F_dada2.clean.fastq.gz"))
filtRs <- file.path(path, "Dada2_Filtered", paste0(sample.names, "_R_dada2.clean.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

filter_output <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
```

Learn errors.
```r
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
```

Sample inference.
```r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

Merge reads.
```r
merged_pairs <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

Get sequence table.
```r
seqtab <- makeSequenceTable(merged_pairs)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)
```

Track reads.
```r
getN <- function(x) sum(getUniques(x))
track <- cbind(filter_output, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(merged_pairs, getN), rowSums(seqtab.nochim), 
               rowSums(seqtab.nochim)/filter_output[,2])
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", 
                     "merged", "nonchim", "surviving reads")
rownames(track) <- sample.names
```

### Assign Taxonomy

Note: A more recent version of the SILVA database has been used (v138.2). The taxonomy assignment was performed on the HPC cluster. The resulting `taxa.rds` object was loaded locally.
```r
red_sea_taxa <- readRDS("./Data/taxa.rds")

red_sea_taxa.print <- red_sea_taxa
rownames(red_sea_taxa.print) <- NULL
head(red_sea_taxa.print)
```

### Extract ASV table
```r
# Give sequence headers manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Write FASTA of final ASV sequences
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "./Data/ASVs.fa")

# Count table: row names without > prefix
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "./Data/ASVs_Counts.tsv", sep="\t", quote=F, col.names=NA)

# Taxonomy table: row names without > prefix
row.names(red_sea_taxa.print) <- sub(">", "", asv_headers)
write.table(red_sea_taxa.print, "./Data/ASVs_Taxonomy.tsv", sep="\t", quote=F, col.names=NA)
```

### Filter mitochondria and chloroplasts

Note: contaminant lists are generated from row names of `red_sea_taxa.print`, which do NOT
have a `>` prefix — consistent with `row.names(asv_tab)`. This ensures the filter matches correctly.
```r
# Filter mitochondria (column 5 = Family) and chloroplasts (column 4 = Order)
mt_ASVs     <- row.names(red_sea_taxa.print)[which(red_sea_taxa.print[, 5] == "Mitochondria")]
chloro_ASVs <- row.names(red_sea_taxa.print)[which(red_sea_taxa.print[, 4] == "Chloroplast")]

# Verify filter works correctly
cat("Mitochondria ASVs found: ", length(mt_ASVs), "\n")
cat("Chloroplast ASVs found:  ", length(chloro_ASVs), "\n")

asv_tab_noMt_noChloro <- asv_tab[!(row.names(asv_tab) %in% c(mt_ASVs, chloro_ASVs)),]

# Write intermediate output files
write.table(asv_tab_noMt_noChloro, "./Data/ASVs_Counts_noMt_noChloro.tsv", 
            sep="\t", quote=F, col.names=NA)

red_sea_taxa.print_noMt_noChloro <- red_sea_taxa.print[!(row.names(red_sea_taxa.print) %in% 
                                                           c(mt_ASVs, chloro_ASVs)),]
write.table(red_sea_taxa.print_noMt_noChloro, "./Data/ASVs_Taxonomy_noMt_noChloro.tsv", 
            sep="\t", quote=F, col.names=NA)
```

### BLAST-validated contaminants

Some taxon assignments cannot be trusted if classification does not go beyond the kingdom level. To remove these, NA entries in the taxon table (= NA at the Kingdom and Phylum level) were checked using BLAST against the NCBI `nt` database to check for eukaryotic contamination. Compared to v138.1, two ASVs (ASV_2167 and ASV_5791) were successfully assigned by SILVA v138.2 and removed from the contaminant list. New contaminants identified in this update include mitochondrial sequences from marine eukaryotes (*Amphimedon queenslandica*, *Callyspongia plicifera*, *Lotharella oceanica*), ASVs with weak BLAST hits (E > 1e-10), and ASVs with no BLAST hits at all.
```r
# Export NA ASVs for BLAST validation
na_asvs <- row.names(red_sea_taxa.print)[is.na(red_sea_taxa.print[, 1]) | 
                                          is.na(red_sea_taxa.print[, 2])]
na_asvs <- na_asvs[na_asvs %in% row.names(asv_tab)]
na_seqs <- asv_seqs[row.names(red_sea_taxa.print) %in% na_asvs]
write(c(rbind(paste0(">", na_asvs), na_seqs)), "./Data/ASVs_NA_forBLAST_v138.2.fa")

# Split into batches of 50 for NCBI BLAST submission
batch_size <- 50
batches <- split(1:length(na_asvs), ceiling(1:length(na_asvs) / batch_size))
for(i in seq_along(batches)){
  idx <- batches[[i]]
  batch_seqs <- c(rbind(paste0(">", na_asvs[idx]), na_seqs[idx]))
  write(batch_seqs, paste0("./Data/BLAST_batch_", i, ".fa"))
}

blast_contaminants <- unique(c(
  # --- Retained from v138.1 (ASV_2167 and ASV_5791 removed: now assigned) ---
  "ASV_12", "ASV_123", "ASV_363", "ASV_499", "ASV_1077",
  "ASV_1703", "ASV_1939", "ASV_2036", "ASV_2209",
  "ASV_2905", "ASV_2955", "ASV_3103", "ASV_3117", "ASV_3238",
  "ASV_4171", "ASV_4296", "ASV_5143", "ASV_5376", "ASV_5382",
  "ASV_5703", "ASV_6449", "ASV_6667", "ASV_6688",
  "ASV_6958", "ASV_6959", "ASV_7036", "ASV_7106", "ASV_7608",
  "ASV_7610", "ASV_7751", "ASV_7753", "ASV_7882", "ASV_8172",
  "ASV_8236", "ASV_8266", "ASV_8830", "ASV_8832", "ASV_9296",
  "ASV_9419", "ASV_9454", "ASV_9460", "ASV_9461", "ASV_9470",
  # --- New v138.2: eukaryotic mitochondrial DNA ---
  # Amphimedon queenslandica mitochondrion (NC_008944.1)
  "ASV_123",
  # Callyspongia plicifera mitochondrion (NC_010206.1)
  "ASV_8333",
  # Lotharella oceanica mitochondrion (NC_029731.1 / NR_041489.2)
  "ASV_7751", "ASV_9167",
  # --- New v138.2: weak E-value BLAST hits (E > 1e-10, likely artefacts) ---
  "ASV_2976", "ASV_3115", "ASV_5109", "ASV_6957",
  "ASV_8881", "ASV_8895", "ASV_8958",
  # --- New v138.2: no BLAST hits (likely artefacts) ---
  "ASV_3032", "ASV_3852", "ASV_3982", "ASV_4450",
  "ASV_5878", "ASV_6763", "ASV_7019", "ASV_7530",
  "ASV_8043", "ASV_8544"
))

# Combine all contaminant ASVs
all_contaminants <- unique(c(mt_ASVs, chloro_ASVs, blast_contaminants))

# Filter tables
asv_tab_clean            <- asv_tab[!(row.names(asv_tab) %in% all_contaminants),]
red_sea_taxa.print_clean <- red_sea_taxa.print[!(row.names(red_sea_taxa.print) %in% 
                                                   all_contaminants),]

# Write output files
write.table(asv_tab_clean,            "./Data/ASVs_Counts_clean_v138.2.tsv", 
            sep="\t", quote=F, col.names=NA)
write.table(red_sea_taxa.print_clean, "./Data/ASVs_Taxonomy_clean_v138.2.tsv", 
            sep="\t", quote=F, col.names=NA)
saveRDS(asv_tab_clean,            "./Data/asv_tab_clean_v138.2.rds")
saveRDS(red_sea_taxa.print_clean, "./Data/taxa_clean_v138.2.rds")
```

### Contaminant summary
```r
contaminant_reads <- sum(asv_tab[row.names(asv_tab) %in% all_contaminants,])
total_reads       <- sum(asv_tab)
cat("=== Contaminant summary ===\n")
cat("Mitochondria removed:       ", length(mt_ASVs),            "ASVs\n")
cat("Chloroplasts removed:       ", length(chloro_ASVs),        "ASVs\n")
cat("BLAST contaminants removed: ", length(blast_contaminants), "ASVs\n")
cat("Total removed:              ", length(all_contaminants),   "ASVs\n")
cat("ASVs before filtering:      ", nrow(asv_tab),              "\n")
cat("ASVs after filtering:       ", nrow(asv_tab_clean),        "\n")
cat("Contaminant reads:          ", contaminant_reads,          "\n")
cat("Total reads:                ", total_reads,                "\n")
cat("Percentage contaminant:     ", round(contaminant_reads / total_reads * 100, 2), "%\n")
cat("Percentage retained:        ", round((1 - contaminant_reads / total_reads) * 100, 2), "%\n")
```

### References

Callahan BJ et al. DADA2: High resolution sample inference from Illumina amplicon data. *Nat Methods* 2016;**13**:581–83. https://doi.org/10.1038/nmeth.3869

Quast C et al. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Res* 2013;**41**:590–96. https://doi.org/10.1093/nar/gks1219