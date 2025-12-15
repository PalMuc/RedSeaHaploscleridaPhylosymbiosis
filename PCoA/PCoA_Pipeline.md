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

### Extract PCoA coordinates
```python
# Get sample scores (coordinates)
pcoa_coords <- data.frame(pcoa_bray$vectors[, 1:3])  # First 3 axes
colnames(pcoa_coords) <- c("PC1", "PC2", "PC3")

# Add metadata
pcoa_coords$Sample <- rownames(pcoa_coords)
pcoa_coords <- merge(pcoa_coords, metadata_clean, by.x = "Sample", by.y = "row.names")

# Save coordinates
write.csv(pcoa_coords, "Beta_diversity_results/pcoa_coordinates.csv", row.names = FALSE)
```

### Basic PCoA plot
```python
plot_pcoa_basic <- ggplot(pcoa_coords, aes(x = PC1, y = PC2, fill = Clade)) +
  geom_point(size = 4, alpha = 0.8, shape = 21, color = "black") +
  scale_fill_manual(values = colour_clades) +
  labs(
    title = "PCoA based on Bray-Curtis distances",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(plot_pcoa_basic)
ggsave("Beta_diversity_results/PCoA_basic.pdf", plot_pcoa_basic, 
       width = 10, height = 7)
```




