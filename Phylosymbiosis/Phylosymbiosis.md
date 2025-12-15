# R-Script Phylosymbiosis Pipeline

This is a modified R-Script following the pipeline as described in Van der Windt et al. 2025.
See original script on: github.com/nielsvanderwindt/2025_Petrosiidae-phylosymbiosis

Before starting this pipeline, I recommend looking at the DADA2 pipeline. In that respective folder, all the files and outputs are further used here.

### Setup
```python
# Load necessary packages
library(phyloseq)
library(microbiome)
library(vegan)
library(ape)
library(phytools)
library(dendextend)
library(phangorn)
library(ggtree)
library(rcompanion)
library(multcompView)
library(picante)
library(pairwiseAdonis)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(DT)
```

### Set working directory
```python
setwd("C:/Path/to/R-Phylosymbiosis/")
```

### Prepare all the data
```python
## Load DADA2 outputs
asv_counts <- read.table("ASVs_Counts_noMt_noChloro.tsv", 
                         sep="\t", header=TRUE, row.names=1)

asv_taxonomy <- read.table("ASVs_Taxonomy_noMt_noChloro.tsv", 
                           sep="\t", header=TRUE, row.names=1)

## Load metadata
sample_metadata <- read.csv("Sample_Information.csv", sep=";")
rownames(sample_metadata) <- sample_metadata$Sample
head(sample_metadata)
rownames(sample_metadata)[1:10] 
all(colnames(asv_counts) %in% rownames(sample_metadata))
all(rownames(sample_metadata) %in% colnames(asv_counts))

## Load phylogenetic tree
asv_tree <- read.tree("ASV_tree_GTR.nwk")

## Create phyloseq object
phyloseq_object <- phyloseq(
  otu_table(asv_counts, taxa_are_rows=TRUE),
  tax_table(as.matrix(asv_taxonomy)),
  sample_data(sample_metadata),
  phy_tree(asv_tree)
)

print(phyloseq_object)

# Save original object
saveRDS(phyloseq_object, "phyloseq_object.rds")
```

