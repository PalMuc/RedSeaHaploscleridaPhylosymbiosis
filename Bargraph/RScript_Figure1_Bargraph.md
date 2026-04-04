# Bargraph Microbiome Red Sea Haplosclerida

Before starting with this pipeline, make sure to have checked the DADA2 pipeline first!

### Setup libraries
```r
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(ape)
library(ggplot2)
```

Make sure the following datasets are loaded:

1. asv_tab_clean (count table)
2. red_sea_taxa.print_clean (taxonomy table)
3. Sample_Information (sample metadata with Clade column)

### Check data
```r
dim(asv_tab_clean)  # ASVs x Samples
dim(red_sea_taxa.print_clean)  # ASVs x Taxonomy ranks
dim(Sample_Information)  # Samples x metadata columns
```

### Create phyloseq object
```r
OTU <- otu_table(asv_tab_clean, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(red_sea_taxa.print_clean))
samples <- sample_data(Sample_Information)
rownames(samples) <- Sample_Information$Sample

physeq <- phyloseq(OTU, TAX, samples)
print(physeq)
```

### Aggregate to class level (not phylum)
This allows us to split Pseudomonadota (formerly Proteobacteria) into Alpha/Gamma/etc.
```r
physeq_class <- tax_glom(physeq, taxrank = "Class")
print(paste("Number of Classes:", ntaxa(physeq_class)))
```

### Transform to relative abundance (%)
```r
physeq_class_rel <- transform_sample_counts(physeq_class, function(x) x/sum(x)*100)
```

### Melt to long format
```r
class_data <- psmelt(physeq_class_rel)
print(paste("Dimensions of class_data:", paste(dim(class_data), collapse = " x ")))
head(class_data)
```

### Check Pseudomonadota classes
```r
pseudomonadota_check <- class_data %>%
  filter(Phylum == "Pseudomonadota") %>%
  group_by(Class) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance))

print("Pseudomonadota Classes:")
print(pseudomonadota_check)
```

### Check Alphaproteobacteria and Gammaproteobacteria
```r
alpha_gamma_check <- class_data %>%
  filter(Class %in% c("Alphaproteobacteria", "Gammaproteobacteria")) %>%
  group_by(Class) %>%
  summarise(
    mean_abund = mean(Abundance),
    n_samples_present = sum(Abundance > 0),
    max_abund = max(Abundance)
  )

print("Alpha vs Gamma presence:")
print(alpha_gamma_check)
```

### Create Phylum_plot column for visualization
```r
# Groups classes into plottable categories
class_data_plot <- class_data %>%
  mutate(Phylum_plot = case_when(
    # Split Pseudomonadota (formerly Proteobacteria) by class
    Class == "Alphaproteobacteria" ~ "Alphaproteobacteria",
    Class == "Gammaproteobacteria" ~ "Gammaproteobacteria",
    Phylum == "Pseudomonadota"     ~ paste0(Class),  # Other Pseudomonadota classes
    # Keep major phyla as-is (SILVA v138.2 names)
    Phylum == "Cyanobacteriota"    ~ "Cyanobacteriota",
    Phylum == "Chloroflexota"      ~ "Chloroflexota",
    Phylum == "Bacteroidota"       ~ "Bacteroidota",
    Phylum == "Acidobacteriota"    ~ "Acidobacteriota",
    Phylum == "Actinomycetota"     ~ "Actinomycetota",
    Phylum == "Planctomycetota"    ~ "Planctomycetota",
    # Everything else
    TRUE ~ "Other"
  ))

table(class_data_plot$Phylum_plot)
```

### Calculate mean abundance per group
```r
plot_summary <- class_data_plot %>%
  group_by(Phylum_plot) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance))

print("Mean abundances for plotting:")
print(plot_summary, n = 20)
```

### Aggregate data per sample
Sum over all classes within each Phylum_plot category.
```r
plot_data <- class_data_plot %>%
  group_by(Sample, Phylum_plot, Clade) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

dim(plot_data)
head(plot_data)
```

### Load phylogeny to get correct sample order
```r
rooted_red_sea_phylo <- read.tree("./Path/to/directory/Rooted_RAxML_RSHaplos_25.tre")

# Get phylogenetic order of samples
phylo_order <- rooted_red_sea_phylo$tip.label

# Factor Sample column in phylogenetic order
plot_data$Sample <- factor(plot_data$Sample, levels = phylo_order)

# Check for NAs (samples not in phylogeny)
sum(is.na(plot_data$Sample))
```

### Define colour palette
```r
phylum_colors <- c(
  "Alphaproteobacteria" = "#5E81AC",    # Blue grey
  "Gammaproteobacteria" = "#BF616A",    # Red
  "Cyanobacteriota"     = "#EBCB8B",    # Gold
  "Bacteroidota"        = "#A3BE8C",    # Green
  "Chloroflexota"       = "#B48EAD",    # Purple
  "Acidobacteriota"     = "#88C0D0",    # Light blue
  "Actinomycetota"      = "#D08770",    # Orange
  "Planctomycetota"     = "#8FBCBB",    # Turquoise
  "Other"               = "#4C566A"     # Dark grey
)
```

### Scale abundance and prepare plot data
```r
plot_data_scaled <- plot_data %>%
  mutate(Abundance_scaled = Abundance / 100)  # Convert percentage to proportion

# Custom legend order (SILVA v138.2 names)
taxa_custom_order <- c(
  "Alphaproteobacteria",
  "Gammaproteobacteria",
  "Cyanobacteriota",
  "Bacteroidota",
  "Chloroflexota",
  "Acidobacteriota",
  "Actinomycetota",
  "Planctomycetota",
  "Other"
)

# Determine major taxa (threshold > 0.25% mean abundance)
major_taxa_only <- plot_summary %>%
  filter(mean_abundance > 0.25, Phylum_plot != "Other") %>%
  pull(Phylum_plot)

print("Taxa to show in plot:")
print(major_taxa_only)

# Reclassify minor taxa as "Other" instead of removing them
# This ensures bars sum to 1.0 for all samples
plot_data_filtered <- plot_data_scaled %>%
  mutate(Phylum_plot = ifelse(Phylum_plot %in% major_taxa_only,
                               as.character(Phylum_plot), "Other")) %>%
  filter(!is.na(Clade)) %>%
  group_by(Sample, Phylum_plot, Clade) %>%
  summarise(Abundance_scaled = sum(Abundance_scaled), .groups = "drop")

# Factor with custom order
plot_data_filtered$Phylum_plot <- factor(plot_data_filtered$Phylum_plot,
                                          levels = taxa_custom_order)

# Check - should be close to 1.0 for all samples
sample_totals <- plot_data_filtered %>%
  group_by(Sample) %>%
  summarise(total = sum(Abundance_scaled))

print("Sample abundance range after filtering:")
print(paste("Min:", round(min(sample_totals$total), 3),
            "Max:", round(max(sample_totals$total), 3)))

# Order clades according to phylogeny
clade_phylo_order <- Sample_Information %>%
  mutate(phylo_position = match(Sample, phylo_order)) %>%
  filter(!is.na(Clade), !is.na(phylo_position)) %>%
  group_by(Clade) %>%
  summarise(first_appearance = min(phylo_position)) %>%
  arrange(first_appearance) %>%
  pull(Clade)

print("Clade order according to phylogeny:")
print(clade_phylo_order)

# Make Clade a factor in phylogenetic order
plot_data_filtered$Clade <- factor(plot_data_filtered$Clade,
                                    levels = clade_phylo_order)
```

### Create the bargraph
```r
Barplot_RS <- ggplot(plot_data_filtered,
                     aes(x = Sample, y = Abundance_scaled, fill = Phylum_plot)) +
  geom_bar(stat = "identity",
           width = 1.0,
           color = "white",
           linewidth = 0.1) +
  # Facet per clade - labels on top
  facet_grid(. ~ Clade,
             scales = "free_x",
             space = "free_x") +
  # Colours (full palette including Other)
  scale_fill_manual(values = phylum_colors,
                    name = "Taxon") +
  # Y-axis
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  # Theme
  theme_minimal() +
  theme(
    # Axes
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    # Facet labels (clade names)
    strip.background = element_rect(fill = "grey", color = NA),
    strip.text = element_text(size = 6, face = "bold", color = "black",
                              margin = margin(3, 0, 3, 0)),
    # Legend
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.5, "cm"),
    # Panel
    panel.spacing    = unit(0.5, "lines"),
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "grey92", color = NA),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  labs(
    x = "Specimens (grouped by host clade)",
    y = "Relative abundance"
  )

print(Barplot_RS)
```

### Save the plot
```r
ggsave("Barplot_RS_v138.2.pdf", plot = Barplot_RS,
       width = 14, height = 5, dpi = 300)
ggsave("Barplot_RS_v138.2.png", plot = Barplot_RS,
       width = 14, height = 5, dpi = 300)
```

### References
Callahan BJ et al. DADA2: High resolution sample inference from Illumina amplicon data. *Nat Methods* 2016;**13**:581–83. https://doi.org/10.1038/nmeth.3869

Quast C et al. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Res* 2013;**41**:590–96. https://doi.org/10.1093/nar/gks1219