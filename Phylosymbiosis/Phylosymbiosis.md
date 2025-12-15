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
# Load DADA2 outputs
asv_counts <- read.table("ASVs_Counts_noMt_noChloro.tsv", 
                         sep="\t", header=TRUE, row.names=1)

asv_taxonomy <- read.table("ASVs_Taxonomy_noMt_noChloro.tsv", 
                           sep="\t", header=TRUE, row.names=1)

# Load metadata
sample_metadata <- read.csv("Sample_Information.csv", sep=";")
rownames(sample_metadata) <- sample_metadata$Sample
head(sample_metadata)
rownames(sample_metadata)[1:10] 
all(colnames(asv_counts) %in% rownames(sample_metadata))
all(rownames(sample_metadata) %in% colnames(asv_counts))

# Load phylogenetic tree
asv_tree <- read.tree("ASV_tree_GTR.nwk")

# Create phyloseq object
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

### Quality control + filtering
```python
# Rarefaction curves
rarecurve(as.data.frame(t(otu_table(phyloseq_object))), 
          step = 100, col = "blue", label = FALSE,
          main = "Rarefaction Curves")

## Check read distribution
print(summary(sample_sums(phyloseq_object)))

hist(sample_sums(phyloseq_object), breaks = 30, col = "steelblue",
     main = "Read counts per sample", 
     xlab = "Number of reads")

## Identify and remove low-quality samples
low_threshold <- 5000
low_samples <- names(sample_sums(phyloseq_object)[sample_sums(phyloseq_object) < low_threshold])
print(low_samples)

# Create filtered dataset
phyloseq_filtered <- prune_samples(sample_sums(phyloseq_object) >= low_threshold, phyloseq_object)
phyloseq_filtered <- prune_taxa(taxa_sums(phyloseq_filtered) > 0, phyloseq_filtered)

print(phyloseq_filtered)
print(summary(sample_sums(phyloseq_filtered)))
```

### Create data subsets
```python
# Rarefied dataset (for alpha diversity)
phyloseq_rarefied <- rarefy_even_depth(phyloseq_filtered, rngseed = 123)
print(phyloseq_rarefied)

# Compositional dataset (for beta diversity)
phyloseq_compositional <- microbiome::transform(phyloseq_filtered, "compositional")
print(phyloseq_compositional)

# Remove samples without clade assignment (for phylosymbiosis)
phyloseq_clades <- subset_samples(phyloseq_filtered, !is.na(Clade))

# Merge by clade
phyloseq_merged <- merge_samples(phyloseq_clades, "Clade")

# Fix sample_data
sample_data(phyloseq_merged)$Clade <- sample_names(phyloseq_merged)
print(phyloseq_merged)
```

### Taxonomic summary
```python
# Phylum-level aggregation
phyloseq_phylum <- tax_glom(phyloseq_filtered, "Phylum")
phylum_abundance <- sort(colSums(t(otu_table(phyloseq_phylum))), decreasing = TRUE)
print(head(phylum_abundance, 15))

# Create taxonomy dataframe with abundance
tax_data <- data.frame(
  as(tax_table(phyloseq_filtered), "matrix"),
  Abundance = taxa_sums(phyloseq_filtered),
  ASV = taxa_names(phyloseq_filtered)
)
```

### Exploratory plots
```python
# ASV abundance distribution
phyloseq_plot1 <- ggplot(tax_data, aes(Abundance)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8) +
  scale_x_log10(labels = scales::comma) +
  labs(title = "ASV Abundance Distribution",
       x = "Abundance (log10 scale)",
       y = "Number of ASVs") +
  theme_bw(base_size = 12)

print(phyloseq_plot1)

# Sample read distribution
sample_data_plot <- data.frame(
  Sample = sample_names(phyloseq_filtered),
  Reads = sample_sums(phyloseq_filtered),
  Clade = sample_data(phyloseq_filtered)$Clade
)

phyloseq_plot2 <- ggplot(sample_data_plot, aes(Reads)) +
  geom_histogram(bins = 30, fill = "coral", alpha = 0.8) +
  labs(title = "Read Distribution per Sample",
       x = "Number of reads",
       y = "Number of samples") +
  theme_bw(base_size = 12)

print(phyloseq_plot2)

# Reads per clade (boxplot)
phyloseq_plot3 <- ggplot(sample_data_plot, aes(x = Clade, y = Reads, fill = Clade)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Read Distribution by Clade",
       x = "Clade",
       y = "Number of reads") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(phyloseq_plot3)
```

### Colour palettes
```python
colour_clades <- c(
  "G01" = "blue3", "G02" = "cadetblue", "G03" = "chocolate4",
  "G04" = "purple3", "G05" = "darkolivegreen", "G06" = "greenyellow",
  "G07" = "cornflowerblue", "G08" = "plum", "G09" = "darkorange1",
  "G10" = "magenta", "G11" = "gold1", "G12" = "cyan",
  "G13" = "bisque2", "G14" = "seagreen1", "G15" = "midnightblue",
  "G16" = "lightcyan3", "G17" = "firebrick2", "G18" = "tan1",
  "G19" = "darkseagreen", "G20" = "maroon2", "G21" = "black",
  "G22" = "red4"
)

# Generate Order-level colours
tax_data$Order_clean <- ifelse(is.na(tax_data$Order) | tax_data$Order == "NA",
                               "Unclassified", tax_data$Order)

tax_data$Full_Order <- paste(
  ifelse(is.na(tax_data$Kingdom), "Unknown", tax_data$Kingdom),
  ifelse(is.na(tax_data$Phylum), "Unclassified", tax_data$Phylum),
  ifelse(is.na(tax_data$Class), "Unclassified", tax_data$Class),
  tax_data$Order_clean,
  sep = " | "
)

unique_orders <- unique(tax_data$Full_Order)
n_orders <- length(unique_orders)

# Generate colours
if(n_orders <= 74) {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  colour_orders <- setNames(col_vector[1:n_orders], unique_orders)
} else {
  colour_orders <- setNames(rainbow(n_orders), unique_orders)
}

colour_orders <- c(colour_orders, "Other" = "grey80")

# Export
write.csv(data.frame(Order = names(colour_orders), Colour = colour_orders),
          "colour_orders.csv", row.names = FALSE)
```

## Beta-diversity

### NMDS plots
```python
# Create output directory
dir.create("Beta_diversity_results", showWarnings = FALSE)

# Calculate NMDS
nmds_bray <- ordinate(phyloseq_compositional, "NMDS", "bray", trymax = 100)

# Save NMDS results
write.csv(nmds_bray$points, "Beta_diversity_results/nmds_bray_points.csv")
saveRDS(nmds_bray, "Beta_diversity_results/nmds_bray.RDS")

# Extract ordination data
nmds_data <- plot_ordination(phyloseq_compositional, nmds_bray, 
                             type = "Samples", justDF = TRUE)

# NMDS plot
plot_nmds <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = Clade), pch = 21, size = 4, alpha = 0.8) +
  scale_fill_manual(values = colour_clades) +
  annotate("text", 
           label = paste0("2D Stress = ", round(nmds_bray$stress, 3)), 
           x = min(nmds_data$NMDS1) + 0.1, 
           y = max(nmds_data$NMDS2) - 0.1,
           size = 5) +
  labs(title = "NMDS based on Bray-Curtis distances",
       fill = "Clade") +
  guides(fill = guide_legend(override.aes = list(size = 7), ncol = 1)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12)
  )

print(plot_nmds)

# Save plot
ggsave(filename = "Plot_NMDS_bray.pdf",
       plot = plot_nmds, 
       device = "pdf", 
       width = 30, height = 20, units = "cm", 
       path = "Beta_diversity_results")

ggsave(filename = "Plot_NMDS_bray.png",
       plot = plot_nmds, 
       device = "png", 
       width = 30, height = 20, units = "cm", dpi = 300,
       path = "Beta_diversity_results")
```

### PERMANOVA
```python
# Get metadata and remove NAs
metadata <- as(sample_data(phyloseq_compositional), "data.frame")
metadata_clean <- metadata[!is.na(metadata$Clade), ]

# Calculate Bray-Curtis distance
samples_with_clades <- rownames(metadata_clean)
phyloseq_comp_clades <- prune_samples(samples_with_clades, phyloseq_compositional)
dist_bray <- distance(phyloseq_comp_clades, method = "bray")

# PERMANOVA test
permanova_result <- adonis2(dist_bray ~ Clade, 
                            data = metadata_clean, 
                            by = NULL,
                            permutations = 999)
print(permanova_result)

# Table
datatable(permanova_result)

# Save results
write.csv(as.data.frame(permanova_result), 
          "Beta_diversity_results/permanova_result_clade.csv")
```


### Pairwise PERMANOVA (post-hoc test)
```python
# Pairwise comparisons between clades
pairwise_result <- pairwise.adonis(
  dist_bray, 
  metadata_clean$Clade,
  p.adjust.m = "BH",
  perm = 999
)
print(head(pairwise_result, 22))

# Interactive table
datatable(pairwise_result)

# Save results
write.csv(as.data.frame(pairwise_result), 
          "Beta_diversity_results/pairwise_permanova_clade.csv")
```

## Mantel Tests

### Data preparation
```python
# Generate microbiome distance matrices

# Use phyloseq_clades (samples with clade info only)
# Bray-Curtis dissimilarity
dist_matrix_bray <- as.matrix(distance(phyloseq_compositional, method = "bray"))

# Weighted UniFrac (phylogenetic + abundance)
dist_matrix_wunifrac <- as.matrix(distance(phyloseq_compositional, method = "wunifrac"))

# Unweighted UniFrac (phylogenetic, presence/absence)
dist_matrix_uunifrac <- as.matrix(distance(phyloseq_compositional, method = "uunifrac"))


# Load host phylogenetic tree and calculate distances

# Load host phylogenetic tree
host_tree <- read.tree("C:/Path/to/directory/Rooted_RAxML_RSHaplos_25.tre")
print(head(host_tree$tip.label, 10))

# Calculate cophenetic distance (patristic distance between tips)
phylo_dist_matrix <- cophenetic(host_tree)

# Save the distance matrix for future use
write.csv(phylo_dist_matrix, 
          "Phylosymbiosis_results/host_phylogenetic_distance_matrix.csv")


# Match sample names between matrices

# Get sample names
micro_samples <- rownames(dist_matrix_bray)
phylo_samples <- rownames(phylo_dist_matrix)

# Find overlapping samples
samples_in_both <- intersect(micro_samples, phylo_samples)

# Check for mismatches
missing_from_phylo <- setdiff(micro_samples, phylo_samples)
missing_from_micro <- setdiff(phylo_samples, micro_samples)

if(length(missing_from_phylo) > 0) {
  cat("\nWARNING: Samples in microbiome but not in phylogeny:\n")
  print(missing_from_phylo)
}

if(length(missing_from_micro) > 0) {
  cat("\nWARNING: Samples in phylogeny but not in microbiome:\n")
  print(missing_from_micro)
}

# Subset matrices to matching samples only
if(length(samples_in_both) < 10) {
  stop("ERROR: Too few matching samples (", length(samples_in_both), 
       "). Check matches between microbiome and phylogeny")
}

# Subset all matrices to matching samples
dist_matrix_bray_match <- dist_matrix_bray[samples_in_both, samples_in_both]
dist_matrix_wunifrac_match <- dist_matrix_wunifrac[samples_in_both, samples_in_both]
dist_matrix_uunifrac_match <- dist_matrix_uunifrac[samples_in_both, samples_in_both]
phylo_dist_matrix_match <- phylo_dist_matrix[samples_in_both, samples_in_both]
```

### Mantel tests
```python
# Mantel test 1: Bray-Curtis vs Host Phylogeny
mantel_bray <- mantel(dist_matrix_bray_match, 
                      phylo_dist_matrix_match, 
                      method = "spearman", 
                      permutations = 9999, 
                      na.rm = TRUE)

# Mantel test 2: Weighted UniFrac vs Host Phylogeny
mantel_wunifrac <- mantel(dist_matrix_wunifrac_match, 
                          phylo_dist_matrix_match, 
                          method = "spearman", 
                          permutations = 9999, 
                          na.rm = TRUE)

# Mantel test 3: Unweighted UniFrac vs Host Phylogeny
mantel_uunifrac <- mantel(dist_matrix_uunifrac_match, 
                          phylo_dist_matrix_match, 
                          method = "spearman", 
                          permutations = 9999, 
                          na.rm = TRUE)

# Create summary table
mantel_results <- data.frame(
  Distance_Method = c("Bray-Curtis", "Weighted UniFrac", "Unweighted UniFrac"),
  Mantel_r = c(mantel_bray$statistic, mantel_wunifrac$statistic, mantel_uunifrac$statistic),
  p_value = c(mantel_bray$signif, mantel_wunifrac$signif, mantel_uunifrac$signif),
  Significance = c(
    ifelse(mantel_bray$signif < 0.001, "***", 
           ifelse(mantel_bray$signif < 0.01, "**",
                  ifelse(mantel_bray$signif < 0.05, "*", "ns"))),
    ifelse(mantel_wunifrac$signif < 0.001, "***",
           ifelse(mantel_wunifrac$signif < 0.01, "**",
                  ifelse(mantel_wunifrac$signif < 0.05, "*", "ns"))),
    ifelse(mantel_uunifrac$signif < 0.001, "***",
           ifelse(mantel_uunifrac$signif < 0.01, "**",
                  ifelse(mantel_uunifrac$signif < 0.05, "*", "ns")))
  )
)

print(mantel_results)

# Save results
write.csv(mantel_results, 
          "Phylosymbiosis_results/mantel_test_results.csv", 
          row.names = FALSE)

# Save matched distance matrices
write.csv(dist_matrix_bray_match, 
          "Phylosymbiosis_results/distance_matrix_bray_matched.csv")
write.csv(dist_matrix_wunifrac_match, 
          "Phylosymbiosis_results/distance_matrix_wunifrac_matched.csv")
write.csv(dist_matrix_uunifrac_match, 
          "Phylosymbiosis_results/distance_matrix_uunifrac_matched.csv")
write.csv(phylo_dist_matrix_match,
          "Phylosymbiosis_results/phylo_distance_matrix_matched.csv")
```



