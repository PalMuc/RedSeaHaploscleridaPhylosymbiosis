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

### PCoA + 95% confidence ellipses
```python
plot_pcoa_ellipse <- ggplot(pcoa_coords, aes(x = PC1, y = PC2, fill = Clade, color = Clade)) +
  stat_ellipse(aes(group = Clade, color = Clade), type = "t", level = 0.95, 
               geom = "polygon", alpha = 0, fill = NA, show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 4, alpha = 0.8, shape = 21, color = "black") +
  scale_fill_manual(values = colour_clades) +
  scale_color_manual(values = colour_clades) +
  labs(
    title = "PCoA with 95% confidence ellipses",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
    fill = "Clade"
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

print(plot_pcoa_ellipse)
ggsave("Beta_diversity_results/PCoA_with_ellipses.pdf", plot_pcoa_ellipse, 
       width = 12, height = 8)
```

### PCoA + centroids
```python
# Calculate centroids per clade
centroids <- pcoa_coords %>%
  group_by(Clade) %>%
  summarise(
    PC1_centroid = mean(PC1),
    PC2_centroid = mean(PC2),
    .groups = "drop"
  )

plot_pcoa_centroids <- ggplot(pcoa_coords, aes(x = PC1, y = PC2, fill = Clade)) +
  geom_point(size = 3, alpha = 0.6, shape = 21, color = "grey30") +
  geom_point(data = centroids, aes(x = PC1_centroid, y = PC2_centroid, fill = Clade),
             size = 6, shape = 23, color = "black", stroke = 1.5) +
  geom_text(data = centroids, aes(x = PC1_centroid, y = PC2_centroid, label = Clade),
            size = 2.5, fontface = "bold") +
  scale_fill_manual(values = colour_clades) +
  labs(
    title = "PCoA with clade centroids",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
    fill = "Clade"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(plot_pcoa_centroids)
ggsave("Beta_diversity_results/PCoA_centroids.pdf", plot_pcoa_centroids, 
       width = 10, height = 8)
```

### Dispersion test
```python
# Test for dispersion among clades
betadisp_result <- betadisper(dist_bray, metadata_clean$Clade)
betadisp_perm <- permutest(betadisp_result, permutations = 999)

print(betadisp_perm)

# Visualise
boxplot(betadisp_result, main = "Distance to centroid by Clade",
        xlab = "Clade", ylab = "Distance to centroid")

# Plot dispersion
plot(betadisp_result, main = "PCoA of Dispersion")
```




