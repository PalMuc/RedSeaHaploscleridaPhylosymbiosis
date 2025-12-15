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
