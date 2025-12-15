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

### Mantel test: Trait vs Microbiome
```python
mantel_trait <- mantel(
  as.dist(micro_dist_all),
  as.dist(trait_dist_matrix),
  method = "spearman",
  permutations = 9999
)

cat("Trait (HMA/LMA) vs Microbiome:\n")
cat("  Mantel r =", round(mantel_trait$statistic, 4), "\n")
cat("  p-value =", mantel_trait$signif, "\n")
if(mantel_trait$signif < 0.001) cat("  ***\n\n")


# Compare: Trait vs Phylogeny as predictors
cat("COMPARISON: Which is a better predictor of microbiome dissimilarity?\n")
cat("  Trait similarity:        r =", round(mantel_trait$statistic, 4), "\n")
cat("  Phylogenetic distance:   r =", round(partial_mantel_phylo$statistic, 4), 
    "(controlling for geography)\n")
cat("  Raw phylogenetic:        r =", round(mantel_bray$statistic, 4), "\n\n")


# Save comparison
trait_vs_phylo <- data.frame(
  Predictor = c("Trait (HMA/LMA)", 
                "Phylogeny (raw)", 
                "Phylogeny (controlling for geography)"),
  Mantel_r = c(mantel_trait$statistic, 
               mantel_bray$statistic, 
               partial_mantel_phylo$statistic),
  p_value = c(mantel_trait$signif, 
              mantel_bray$signif, 
              partial_mantel_phylo$signif)
)

write.csv(trait_vs_phylo,
          "Phylosymbiosis_results/trait_vs_phylogeny_comparison.csv",
          row.names = FALSE)
```

## PERMANOVA

### Trait effect
```python
# Subset to samples with trait classification (exclude > Intermediate)
samples_HMA_LMA <- rownames(metadata_all)[metadata_all$Microbial_Type %in% c("HMA", "LMA")]

phyloseq_HMA_LMA <- prune_samples(samples_HMA_LMA, phyloseq_compositional)

# Get distance for these samples
dist_HMA_LMA <- distance(phyloseq_HMA_LMA, method = "bray")
metadata_HMA_LMA <- as(sample_data(phyloseq_HMA_LMA), "data.frame")

# PERMANOVA: Trait effect
# How many of the total microbiome variation can be explained by HMA vs. LMA classification?
permanova_trait <- adonis2(
  dist_HMA_LMA ~ Microbial_Type,
  data = metadata_HMA_LMA,
  permutations = 999
)

print(permanova_trait)

# Save results
write.csv(as.data.frame(permanova_trait),
          "Phylosymbiosis_results/permanova_trait_HMA_LMA.csv")
```

### Clade effect (for comparison)
```python
# How many of the microbiome variation observed can be explained by the specific clade identity?
permanova_clade_HMA_LMA <- adonis2(
  dist_HMA_LMA ~ Clade,
  data = metadata_HMA_LMA,
  permutations = 999
)

print(permanova_clade_HMA_LMA)

# Save results
write.csv(as.data.frame(permanova_clade_HMA_LMA),
          "Phylosymbiosis_results/permanova_clade_within_HMA_LMA.csv")
```

## Variance Partitioning Analysis

How many of the microbiome variation can be explained by traits, phylogeny and geography?
How much overlap is there between these variables?

### Partitioning microbiome variance among Trait, Phylogeny, and Geography
```python
# Use samples that have all data and are HMA or LMA
samples_varpart <- intersect(samples_all_data, samples_HMA_LMA)

# Subset matrices
micro_dist_vp <- micro_dist_all[samples_varpart, samples_varpart]
phylo_dist_vp <- phylo_dist_all[samples_varpart, samples_varpart]
geo_dist_vp <- geo_dist_all[samples_varpart, samples_varpart]
trait_dist_vp <- trait_dist_matrix[samples_varpart, samples_varpart]

# Convert to dist objects
micro_dist_vp <- as.dist(micro_dist_vp)

# Get metadata for these samples
metadata_vp <- metadata_all[samples_varpart, ]

# Use PCoA coordinates for continuous predictors
phylo_pcoa <- cmdscale(phylo_dist_vp, k = 3)
colnames(phylo_pcoa) <- paste0("Phylo_PC", 1:3)

geo_pcoa <- cmdscale(geo_dist_vp, k = 3)
colnames(geo_pcoa) <- paste0("Geo_PC", 1:3)

# Combine with metadata
metadata_vp_extended <- cbind(metadata_vp, phylo_pcoa, geo_pcoa)

# Sequential R² PERMANOVA for variance partitioning
# Pure effects (single factor models)
trait_only <- adonis2(micro_dist_vp ~ Microbial_Type, 
                      data = metadata_vp_extended, permutations = 999)

phylo_only <- adonis2(micro_dist_vp ~ Phylo_PC1 + Phylo_PC2 + Phylo_PC3,
                      data = metadata_vp_extended, permutations = 999)

geo_only <- adonis2(micro_dist_vp ~ Geo_PC1 + Geo_PC2 + Geo_PC3,
                    data = metadata_vp_extended, permutations = 999)

# Full model
full_model <- adonis2(micro_dist_vp ~ Microbial_Type + 
                        Phylo_PC1 + Phylo_PC2 + Phylo_PC3 + 
                        Geo_PC1 + Geo_PC2 + Geo_PC3,
                      data = metadata_vp_extended, permutations = 999)
```

### Extract R² values properly
```python
# The "Model" row already contains the TOTAL R² for all terms
R2_trait_only <- trait_only$R2[1]  # Row 1 = Model (already total)
R2_phylo_only <- phylo_only$R2[1]  # Row 1 = Model (already total for all 3 PCs!)
R2_geo_only <- geo_only$R2[1]      # Row 1 = Model (already total for all 3 PCs!)
R2_full <- full_model$R2[1]        # Row 1 = Model (already total for all 7 terms!)

# Overlap estimate
overlap_estimate <- (R2_trait_only + R2_phylo_only + R2_geo_only) - R2_full

# Save variance partitioning results
varpart_summary <- data.frame(
  Factor = c("Trait (HMA/LMA)", "Phylogeny (3 PCs)", "Geography (3 PCs)", 
             "Full model", "Shared/Overlap", "Unexplained"),
  R_squared = c(R2_trait_only, R2_phylo_only, R2_geo_only, 
                R2_full, overlap_estimate, 1 - R2_full),
  Percentage = round(c(R2_trait_only, R2_phylo_only, R2_geo_only, 
                       R2_full, overlap_estimate, 1 - R2_full) * 100, 2)
)


write.csv(varpart_summary,
          "Phylosymbiosis_results/variance_partitioning_summary.csv",
          row.names = FALSE)


# Variance Partitioning Summary
print(varpart_summary)
```

## Within vs. between group comparisons

Are microbiomes from samples with the same trait type (HMA vs HMA or LMA vs LMA) more similar to each other, 
than microbiomes from samples that have different trait types (HMA vs LMA)?

### Calculate mean Bray-Curtis distances
```python
# Convert to long format
micro_dist_long <- as.data.frame(as.matrix(micro_dist_all))
micro_dist_long$Sample1 <- rownames(micro_dist_long)
micro_dist_long <- pivot_longer(micro_dist_long, 
                                cols = -Sample1, 
                                names_to = "Sample2", 
                                values_to = "Bray_Curtis")

# Add trait info
micro_dist_long$Trait1 <- metadata_all[micro_dist_long$Sample1, "Microbial_Type"]
micro_dist_long$Trait2 <- metadata_all[micro_dist_long$Sample2, "Microbial_Type"]

# Remove self-comparisons
micro_dist_long <- micro_dist_long[micro_dist_long$Sample1 != micro_dist_long$Sample2, ]

# Classify comparisons
micro_dist_long$Comparison_Type <- case_when(
  micro_dist_long$Trait1 == "HMA" & micro_dist_long$Trait2 == "HMA" ~ "Within HMA",
  micro_dist_long$Trait1 == "LMA" & micro_dist_long$Trait2 == "LMA" ~ "Within LMA",
  (micro_dist_long$Trait1 == "HMA" & micro_dist_long$Trait2 == "LMA") |
    (micro_dist_long$Trait1 == "LMA" & micro_dist_long$Trait2 == "HMA") ~ "Between HMA-LMA",
  TRUE ~ "Other"
)

# Calculate summaries
comparison_summary <- micro_dist_long %>%
  filter(Comparison_Type != "Other") %>%
  group_by(Comparison_Type) %>%
  summarise(
    Mean_Distance = mean(Bray_Curtis, na.rm = TRUE),
    SD_Distance = sd(Bray_Curtis, na.rm = TRUE),
    Median_Distance = median(Bray_Curtis, na.rm = TRUE),
    N = n()
  )

# Mean Bray-Curtis distances
print(comparison_summary)
```

## Kruskal-Wallis + Post-hoc Dunn Tests

Are the means of these groups significantly different, or is this a difference by coincidence?
(Non-parametric equivalent of a one-way ANOVA) Bray-Curtis values are not normally distributed 
Hence, Kruskal-Wallis test is a better option.

```python
kruskal_result <- kruskal.test(Bray_Curtis ~ Comparison_Type, 
                               data = filter(micro_dist_long, Comparison_Type != "Other"))
                               
cat("\nKruskal-Wallis test:\n")
cat("  Chi-squared =", round(kruskal_result$statistic, 3), "\n")
cat("  df =", kruskal_result$parameter, "\n")
cat("  p-value =", format.pval(kruskal_result$p.value, digits = 3), "\n\n")
```

Kruskal-Wallis test tells us that there is a difference, but not where.
For this, a Post-hoc Dunn test is necessary!

```python
if(kruskal_result$p.value < 0.05) {
  library(FSA)
  dunn_result <- dunnTest(Bray_Curtis ~ Comparison_Type,
                          data = filter(micro_dist_long, Comparison_Type != "Other"),
                          method = "bonferroni")
  cat("Post-hoc Dunn test (Bonferroni correction):\n")
  print(dunn_result$res)
}
```


                               
                               






