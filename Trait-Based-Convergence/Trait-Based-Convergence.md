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

