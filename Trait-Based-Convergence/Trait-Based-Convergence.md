# R-Script for Trait-Based vs. Phylogenetic Filtering 

Before starting with this script, please check the DADA2 + Phylosymbiosis + PCoA pipelines in advance.


### Setup libraries
```python
# Load additional packages
library(vegan)
library(geosphere)  
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)  
library(FSA)
```

### Load data
```python
# Load sample data from Excel (same as used for maps)
sample_data_full <- read_excel("C:/Path/to/directory/R-Phylosymbiosis/Supplementary_Table_1.xlsx")

# Rename ID column to Sample for consistency
sample_data_full <- sample_data_full %>%
  rename(Sample = ID)

# Fix coordinate issues (same correction as in mapping script)
sample_data_full <- sample_data_full %>%
  mutate(Longitude = ifelse(Longitude > 100, Longitude / 1000000, Longitude))

# Remove any samples with missing or invalid coordinates
sample_data_full <- sample_data_full %>%
  filter(!is.na(Latitude), !is.na(Longitude),
         Longitude >= 30, Longitude <= 50,
         Latitude >= 10, Latitude <= 32)
```

### Merge with phyloseq metadata
```python
# Get metadata from phyloseq object
metadata_phyloseq <- as(sample_data(phyloseq_compositional), "data.frame")
metadata_phyloseq$Sample <- rownames(metadata_phyloseq)

# Merge coordinates into phyloseq metadata
metadata_analysis <- metadata_phyloseq %>%
  left_join(sample_data_full %>% select(Sample, Longitude, Latitude, Clade),
            by = "Sample",
            suffix = c("", "_excel"))

# Use Clade from excel if missing in phyloseq
if("Clade_excel" %in% colnames(metadata_analysis)) {
  metadata_analysis <- metadata_analysis %>%
    mutate(Clade = coalesce(Clade, Clade_excel)) %>%
    select(-Clade_excel)
}

rownames(metadata_analysis) <- metadata_analysis$Sample
```

### Prepare HMA/LMA classification
```python
# Based on microbiome composition patterns from manuscript:
# HMA clades: G01, G06, G12 (Chloroflexi-Acidobacteria enriched)
# LMA clades: G03, G04, G05, G07, G09, G15, G16 (Cyanobacteria dominated)
# Intermediate/Other: remaining clades

# Add HMA/LMA classification
metadata_analysis$Microbial_Type <- case_when(
  metadata_analysis$Clade %in% c("G01", "G06", "G12") ~ "HMA",
  metadata_analysis$Clade %in% c("G03", "G04", "G05", "G07", "G09", "G15", "G16") ~ "LMA",
  TRUE ~ "Intermediate"
)

# Check distribution
print(table(metadata_analysis$Microbial_Type))

# Add back to phyloseq object
sample_data(phyloseq_compositional)$Microbial_Type <- metadata_analysis$Microbial_Type

# Save classification
write.csv(metadata_analysis[, c("Sample", "Clade", "Microbial_Type", "Longitude", "Latitude")],
          "Phylosymbiosis_results/HMA_LMA_classification_with_coords.csv",
          row.names = FALSE)
```

### Calculate geographic distances between all samples
```python
# Extract coordinates - now from merged data
coords <- metadata_analysis[, c("Longitude", "Latitude")]
rownames(coords) <- metadata_analysis$Sample

# Remove samples without coordinates
coords_complete <- coords[complete.cases(coords), ]

# Calculate geographic distance matrix (in km)
geo_dist_matrix <- distm(coords_complete, fun = distHaversine) / 1000  # Convert to km
rownames(geo_dist_matrix) <- rownames(coords_complete)
colnames(geo_dist_matrix) <- rownames(coords_complete)

# Save geographic distances
write.csv(geo_dist_matrix,
          "Phylosymbiosis_results/geographic_distance_matrix_km.csv")
```

### Match all matrices (i.e., phylogeny, microbiome, geography)
```python
# Find samples present in all three datasets
samples_with_phylo <- rownames(phylo_dist_matrix)
samples_with_micro <- rownames(dist_matrix_bray)
samples_with_geo <- rownames(geo_dist_matrix)

# Intersection
samples_all_data <- Reduce(intersect, list(samples_with_phylo, 
                                           samples_with_micro, 
                                           samples_with_geo))

if(length(samples_all_data) < 50) {
  warning("Only ", length(samples_all_data), " samples with complete data. Consider checking data completeness.")
}

# Subset all matrices to common samples
phylo_dist_all <- phylo_dist_matrix[samples_all_data, samples_all_data]
micro_dist_all <- dist_matrix_bray[samples_all_data, samples_all_data]
geo_dist_all <- geo_dist_matrix[samples_all_data, samples_all_data]

# Get matches metadata
metadata_all <- metadata_analysis[samples_all_data, ]
```

## Partial Mantel Tests

### Testing phylogeny vs geography effects on microbiome
```python
# Test 1: What is the effect of phylogeny on the microbiome, controlling for geography
partial_mantel_phylo <- mantel.partial(
  as.dist(micro_dist_all),
  as.dist(phylo_dist_all),
  as.dist(geo_dist_all),
  method = "spearman",
  permutations = 9999
)


# Test 2: What is the effect of geography on the microbiome, controlling for phylogeny
partial_mantel_geo <- mantel.partial(
  as.dist(micro_dist_all),
  as.dist(geo_dist_all),
  as.dist(phylo_dist_all),
  method = "spearman",
  permutations = 9999
)


# Test 3: Simple Mantel for Geography (for comparison)
# Is there a correlation between geography and the microbiome?
simple_mantel_geo <- mantel(
  as.dist(micro_dist_all),
  as.dist(geo_dist_all),
  method = "spearman",
  permutations = 9999
)


# Save results
partial_mantel_results <- data.frame(
  Test = c("Phylogeny | Geography", 
           "Geography | Phylogeny",
           "Geography (simple)",
           "Phylogeny (simple, from earlier)"),
  Mantel_r = c(partial_mantel_phylo$statistic,
               partial_mantel_geo$statistic,
               simple_mantel_geo$statistic,
               mantel_bray$statistic),
  p_value = c(partial_mantel_phylo$signif,
              partial_mantel_geo$signif,
              simple_mantel_geo$signif,
              mantel_bray$signif),
  Interpretation = c(
    "Phylogenetic effect after removing geography",
    "Geographic effect after removing phylogeny",
    "Raw geographic effect",
    "Raw phylogenetic effect"
  )
)

write.csv(partial_mantel_results,
          "Phylosymbiosis_results/partial_mantel_results.csv",
          row.names = FALSE)

print(partial_mantel_results[, 1:3])
```

### Trait (HMA vs. LMA) vs. phylogeny comparison
```python
# Create trait distance matrix
# 0 if same type (HMA-HMA or LMA-LMA), 1 if different type
trait_vector <- metadata_all$Microbial_Type
names(trait_vector) <- rownames(metadata_all)

# Create distance matrix
n_samples <- length(trait_vector)
trait_dist_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
rownames(trait_dist_matrix) <- names(trait_vector)
colnames(trait_dist_matrix) <- names(trait_vector)

for(i in 1:n_samples) {
  for(j in 1:n_samples) {
    if(trait_vector[i] != trait_vector[j]) {
      trait_dist_matrix[i, j] <- 1
    }
  }
}

# Save trait distance matrix
write.csv(trait_dist_matrix,
          "Phylosymbiosis_results/trait_distance_matrix.csv")
```






