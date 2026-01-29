# R-Script Tanglegram

Before starting this pipeline, make sure to check the DADA2 and Phylosymbiosis pipelines first!


### Prepare ultrametric tree of the host
```python
if(!is.ultrametric(host_tree_pruned)) {
  cat("Making tree ultrametric...\n")
  host_tree_ultra <- force.ultrametric(host_tree_pruned, method = "nnls")
} else {
  host_tree_ultra <- host_tree_pruned
}

# Quick check plot
pdf("Phylosymbiosis_results/Host_tree_ultrametric_check.pdf", width = 8, height = 10)
plot(host_tree_ultra, main = "Host Phylogenetic Tree (Ultrametric)", cex = 0.6)
dev.off()
```

### Check label compatibility
```python
# Get labels from host tree
host_labels <- host_tree_ultra$tip.label
# Get labels from microbiome dendrograms
micro_labels_bray <- micro_dend_bray$labels

# Find common samples
common_samples <- intersect(host_labels, micro_labels_bray)
# Check for mismatches
if(length(common_samples) < length(host_labels)) {
  cat("\nWARNING: Some samples don't match. Pruning to common samples.\n")
  missing_from_micro <- setdiff(host_labels, micro_labels_bray)
  missing_from_host <- setdiff(micro_labels_bray, host_labels)
  
  if(length(missing_from_micro) > 0) {
    cat("In host but not microbiome:", length(missing_from_micro), "\n")
  }
  if(length(missing_from_host) > 0) {
    cat("In microbiome but not host:", length(missing_from_host), "\n")
  }
}
```

### Rebuild dendrograms with matching samples
```python
# Prune host tree to common samples
host_tree_matched <- keep.tip(host_tree_ultra, common_samples)

# Rebuild microbiome dendrograms from distance matrices with common samples
dist_bray_common <- as.dist(dist_matrix_bray_match[common_samples, common_samples])
dist_wunifrac_common <- as.dist(dist_matrix_wunifrac_match[common_samples, common_samples])
dist_uunifrac_common <- as.dist(dist_matrix_uunifrac_match[common_samples, common_samples])

# Create new dendrograms
micro_dend_bray_matched <- hclust(dist_bray_common, method = "average")
micro_dend_wunifrac_matched <- hclust(dist_wunifrac_common, method = "average")
micro_dend_uunifrac_matched <- hclust(dist_uunifrac_common, method = "average")

## Add clade info to labels
# Get metadata for matching samples
metadata_match <- as(sample_data(phyloseq_compositional), "data.frame")
metadata_match <- metadata_match[common_samples, ]

# Add clade info to host tree labels
host_tree_matched$tip.label <- paste0(
  host_tree_matched$tip.label, 
  " (", 
  metadata_match[host_tree_matched$tip.label, "Clade"],
  ")"
)

# Add same clade info to microbiome dendrogram labels
micro_dend_bray_matched$labels <- paste0(
  micro_dend_bray_matched$labels, 
  " (", 
  metadata_match[micro_dend_bray_matched$labels, "Clade"],
  ")"
)

micro_dend_wunifrac_matched$labels <- paste0(
  micro_dend_wunifrac_matched$labels, 
  " (", 
  metadata_match[micro_dend_wunifrac_matched$labels, "Clade"],
  ")"
)

micro_dend_uunifrac_matched$labels <- paste0(
  micro_dend_uunifrac_matched$labels, 
  " (", 
  metadata_match[micro_dend_uunifrac_matched$labels, "Clade"],
  ")"
)
```

### Convert to dendrograms and create dendlists
```python
# Convert host tree to dendrogram
host_dend_matched <- as.dendrogram(as.hclust(host_tree_matched))

# Convert microbiome hclust to dendrograms
micro_dend_bray_dend <- as.dendrogram(micro_dend_bray_matched)
micro_dend_wunifrac_dend <- as.dendrogram(micro_dend_wunifrac_matched)
micro_dend_uunifrac_dend <- as.dendrogram(micro_dend_uunifrac_matched)
```

### Create dendlists and untangle
```python
# Create dendlists
dends_bray <- dendlist(host_dend_matched, micro_dend_bray_dend)
dends_bray_untangled <- untangle(dends_bray, method = "step2side")

dends_wunifrac <- dendlist(host_dend_matched, micro_dend_wunifrac_dend)
dends_wunifrac_untangled <- untangle(dends_wunifrac, method = "step2side")

dends_uunifrac <- dendlist(host_dend_matched, micro_dend_uunifrac_dend)
dends_uunifrac_untangled <- untangle(dends_uunifrac, method = "step2side")

# Calculate entanglement
entanglement_bray_final <- entanglement(dends_bray_untangled)
entanglement_wunifrac_final <- entanglement(dends_wunifrac_untangled)
entanglement_uunifrac_final <- entanglement(dends_uunifrac_untangled)

cat("Bray-Curtis:", round(entanglement_bray_final, 4), "\n")
cat("Weighted UniFrac:", round(entanglement_wunifrac_final, 4), "\n")
cat("Unweighted UniFrac:", round(entanglement_uunifrac_final, 4), "\n")
```

### Plot tanglegram 
```python
tanglegram(
  dends_bray_untangled, 
  highlight_distinct_edges = FALSE, 
  common_subtrees_color_lines = FALSE, 
  highlight_branches_lwd = FALSE, 
  lwd = 1,
  main_left = "Host Phylogeny",
  main_right = "Microbiome Dendrogram (Bray-Curtis)",
  main = paste0("Entanglement = ", round(entanglement_bray_final, 3)),
  cex_main = 1,
  columns_width = c(2, 2, 2),
  lab.cex = 0.45,
  margin_inner = 3,
  margin_outer = 1
)
```

### References

Van der Windt N et al. Host evolutionary history drives prokaryotic diversity in the globally distributed sponge family Petrosiidae. Mol Ecol 2025;34:e70186. https://doi.org/10.1111/mec.70186
