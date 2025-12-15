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
