# R-Script PCoA Pipeline

For this script you'll need some of the files or variables generated in the Phylosymbiosis pipeline.
Please look at the respective folder for more information.

### Setup
```python
# Install the necessary libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
```

### Prepare data
```python
# Get metadata and remove NA's (use same clean dataset as PERMANOVA)
metadata <- as(sample_data(phyloseq_compositional), "data.frame")
metadata_clean <- metadata[!is.na(metadata$Clade), ]

# Subset phyloseq object to samples with clades
samples_with_clades <- rownames(metadata_clean)
phyloseq_comp_clades <- prune_samples(samples_with_clades, phyloseq_compositional)

# Calculate Bray-Curtis distance matrix
dist_bray <- distance(phyloseq_comp_clades, method = "bray")
```

### Run PCoA
```python
# Perform PCoA ordination
pcoa_bray <- ordinate(phyloseq_comp_clades, method = "PCoA", distance = "bray")

# Extract eigenvalues and calculate % variance explained
eigenvalues <- pcoa_bray$values$Eigenvalues
variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Print variance explained by first 5 axes
for(i in 1:5) {
  cat(paste0("PC", i, ": ", round(variance_explained[i], 2), "%\n"))
}

# Calculate cumulative variance
cumulative_var <- cumsum(variance_explained)
cat(paste0("\nCumulative variance PC1-2: ", round(cumulative_var[2], 2), "%\n"))
cat(paste0("Cumulative variance PC1-3: ", round(cumulative_var[3], 2), "%\n"))
```





