# R-Script for Trait-Based vs. Phylogenetic Filtering Prokaryotes in Sponges

Before starting with this script, please check the DADA2 + Phylosymbiosis pipelines in advance.


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
# HMA clades: G01, G06, G12, G22 (Chloroflexi-Acidobacteria enriched)
# LMA clades: G03, G04, G05, G07, G08, G09, G10, G11 (Cyanobacteria dominated)
# Intermediate/Other: remaining clades

# Add HMA/LMA classification
metadata_analysis$Microbial_Type <- case_when(
  metadata_analysis$Clade %in% c("G01", "G06", "G12", "G22") ~ "HMA",
  metadata_analysis$Clade %in% c("G03", "G04", "G05", "G07", "G08", "G09", "G10", "G11") ~ "LMA",
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

