# Configuration Heatmap Red Sea Haplosclerida-ASVs 

Before starting with this pipeline, please check the Phylosymbiosis pipeline first!

### Setup libraries
```python
library(matrixStats)
library(circlize)
library(ComplexHeatmap)
```

### Create data matrix for heatmap
```python
# Concert to matrix and calculate total abundance ASVs
asv_mat <- as.matrix(asv_tab_noMt_noChloro)

# Calculate total abundance each ASV
asv_totals <- rowSums(asv_mat)

# Select those 100 ASVs with highest total abundance
top_asvs <- names(sort(asv_totals, decreasing = TRUE)[1:100])

# Subset matrix
asv_top <- asv_mat[top_asvs, ]
asv_rel <- prop.table(asv_top, 2)     # normalise for each column
asv_rel_sqrt <- sqrt(asv_rel)         # transform to square-root
asv_heat <- asv_rel_sqrt
```

### Transpose matrix
```python
# Check structure before transpose
print(paste("prior transpose Rows:", nrow(asv_heat), "Columns:", ncol(asv_heat)))
print(head(rownames(asv_heat)))
print(head(colnames(asv_heat)))

# Transpose Sponges as rows, ASVs as columns
asv_heat <- t(asv_heat)

# Check structure after transpose
print(paste("after transpose Rows:", nrow(asv_heat), "Columns:", ncol(asv_heat)))
print(head(rownames(asv_heat)))
print(head(colnames(asv_heat)))

# Order samples (rows) according to phylogenetic tree
sample_order <- rooted_red_sea_phylo$tip.label
sample_order <- sample_order[sample_order %in% rownames(asv_heat)]  # Samples in rows
asv_heat <- asv_heat[sample_order, ]

# Order ASVs (columns) according to excel file
asv_order <- ASV_Order_Heatmap$ASV  # Use as ASV column
asv_order <- asv_order[asv_order %in% colnames(asv_heat)] # Filter for only those ASVs present

# Check number of matches
print(paste(length(asv_order)))
print(paste(ncol(asv_heat)))

# Sort the columns according to ASV order
asv_heat <- asv_heat[, asv_order]

# Verify
print(head(colnames(asv_heat), 10))
```

### Create Heatmap
```python
heat_colors <- colorRamp2(
  c(0, 0.001, max(asv_heat)),           # Breakpoints: 0, min:>0, max:1
  c("#FFFFFF", "#9ECAE1", "#6A1B9A")    # Colours: white, lightblue, deep purple
)

Heatmap(
  asv_heat,
  name = "Rel. abundance",
  col = heat_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_side = "top",           # ASVs on top
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(title = "square root of rel. abundance (%)")
)
```










