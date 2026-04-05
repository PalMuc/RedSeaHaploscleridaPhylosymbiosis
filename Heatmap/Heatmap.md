# Configuration Heatmap Red Sea Haplosclerida-ASVs

Before starting with this pipeline, please check the Phylosymbiosis pipeline first!

### Setup libraries
```r
library(matrixStats)
library(circlize)
library(ComplexHeatmap)
```

### Create data matrix for heatmap
```r
# Convert to matrix and calculate total abundance ASVs
asv_mat <- as.matrix(asv_tab_clean)

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

### Generate ASV order table
This table is used to order the ASV columns in the heatmap by taxonomy (Phylum > Class > Order).
It is generated fresh from the current top 100 ASVs and `red_sea_taxa.print_clean` to ensure
consistency with the filtered dataset.
```r
asv_order_new <- data.frame(
  ASV    = top_asvs,
  Phylum = red_sea_taxa.print_clean[top_asvs, 2],
  Class  = red_sea_taxa.print_clean[top_asvs, 3],
  Order  = red_sea_taxa.print_clean[top_asvs, 4]
)

# Sort by Phylum and Class for logical grouping in heatmap
asv_order_new <- asv_order_new[order(asv_order_new$Phylum,
                                      asv_order_new$Class,
                                      asv_order_new$Order), ]

write.csv(asv_order_new, "ASV_Order_Heatmap_v138.2.csv", row.names = FALSE)
cat("Saved ASV_Order_Heatmap_v138.2.csv with", nrow(asv_order_new), "ASVs\n")
```

### Transpose matrix
```r
# Check structure before transpose
print(paste("prior transpose Rows:", nrow(asv_heat), "Columns:", ncol(asv_heat)))
print(head(rownames(asv_heat)))
print(head(colnames(asv_heat)))

# Transpose: Sponges as rows, ASVs as columns
asv_heat <- t(asv_heat)

# Check structure after transpose
print(paste("after transpose Rows:", nrow(asv_heat), "Columns:", ncol(asv_heat)))
print(head(rownames(asv_heat)))
print(head(colnames(asv_heat)))

# Order samples (rows) according to phylogenetic tree
sample_order <- rooted_red_sea_phylo$tip.label
sample_order <- sample_order[sample_order %in% rownames(asv_heat)]
asv_heat <- asv_heat[sample_order, ]

# Load ASV order table (generated from red_sea_taxa.print_clean, SILVA v138.2)
ASV_Order_Heatmap <- read.csv("ASV_Order_Heatmap_v138.2.csv")

# Order ASVs (columns) according to taxonomy table
asv_order <- ASV_Order_Heatmap$ASV
asv_order <- asv_order[asv_order %in% colnames(asv_heat)]

# Check number of matches
print(paste(length(asv_order)))
print(paste(ncol(asv_heat)))

# Sort the columns according to ASV order
asv_heat <- asv_heat[, asv_order]

# Verify
print(head(colnames(asv_heat), 10))
```

### Create Heatmap
```r
heat_colors <- colorRamp2(
  c(0, 0.001, max(asv_heat)),           # Breakpoints: 0, min:>0, max:1
  c("#FFFFFF", "#9ECAE1", "#6A1B9A")    # Colours: white, lightblue, deep purple
)

heatmap_plot <- Heatmap(
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

draw(heatmap_plot)
```

### Save Heatmap
```r
pdf("Heatmap_RS_v138.2.pdf", width = 14, height = 8)
draw(heatmap_plot)
dev.off()

svg("Heatmap_RS_v138.2.svg", width = 14, height = 8)
draw(heatmap_plot)
dev.off()
```