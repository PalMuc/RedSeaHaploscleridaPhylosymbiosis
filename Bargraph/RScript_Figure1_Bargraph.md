# Bargraph Microbiome Red Sea Haplosclerida

Before starting with this pipeline, make sure to have checked the DADA2 pipeline first!

### Setup libraries
```python
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(ape)
library(ggplot2)
```
Make sure the following datasets are loaded:
1. asv_tab_noMt_noChloro (count table)
2. red_sea_taxa.print_noMt_noChloro (taxonomy table)
3. sample_info (sample metadata with Clade column)

### Check data
```python
dim(asv_tab_noMt_noChloro)  #ASVs x Samples
dim(red_sea_taxa.print_noMt_noChloro)  #ASVs x Taxonomy ranks
dim(sample_info)  # Samples x metadata columns
```

### Create phyloseq object
```python
OTU <- otu_table(asv_tab_noMt_noChloro, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(red_sea_taxa.print_noMt_noChloro))
samples <- sample_data(sample_info)
rownames(samples) <- sample_info$Sample

physeq <- phyloseq(OTU, TAX, samples)
print(physeq)
```

### Aggregate to CLASS level (not Phylum)
This allows us to split Proteobacteria into Alpha/Gamma/etc.
```python
physeq_class <- tax_glom(physeq, taxrank = "Class")
print(paste("Number of Classes:", ntaxa(physeq_class)))
```

### Transform to relative abundance (%)
```python
physeq_class_rel <- transform_sample_counts(physeq_class, function(x) x/sum(x)*100)
```

### Melt to long format
```python
class_data <- psmelt(physeq_class_rel)
print(paste("Dimensions of class_data:", paste(dim(class_data), collapse = " x ")))
head(class_data)
```

### Check what we have for Proteobacteria
```python
proteo_check <- class_data %>%
  filter(Phylum == "Proteobacteria") %>%
  group_by(Class) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance))

print("Proteobacteria Classes:")
print(proteo_check)
```

### Create Phylum_plot column for visualization
```python
# This groups Classes into plottable categories
class_data_plot <- class_data %>%
  mutate(Phylum_plot = case_when(
    # Split Proteobacteria by Class
    Class == "Alphaproteobacteria" ~ "Alphaproteobacteria",
    Class == "Gammaproteobacteria" ~ "Gammaproteobacteria",
    Class == "Deltaproteobacteria" ~ "Deltaproteobacteria",
    Phylum == "Proteobacteria" ~ paste0(Class),  # Other Proteobacteria classes
    # Keep major phyla as-is
    Phylum == "Cyanobacteria" ~ "Cyanobacteria",
    Phylum == "Chloroflexi" ~ "Chloroflexi",
    Phylum == "Bacteroidota" ~ "Bacteroidota",
    Phylum == "Acidobacteriota" ~ "Acidobacteriota",
    Phylum == "Actinobacteriota" ~ "Actinobacteriota",
    Phylum == "Planctomycetota" ~ "Planctomycetota",
    # Everything else = Other
    TRUE ~ "Other"
  ))

table(class_data_plot$Phylum_plot)

```

### Calculate mean abundance per group
```python
plot_summary <- class_data_plot %>%
  group_by(Phylum_plot) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance))

print("Mean abundances for plotting:")
print(plot_summary, n = 20)
```

### Check if we have Alpha + Gamma
```python
alpha_gamma_check <- class_data_plot %>%
  filter(Phylum_plot %in% c("Alphaproteobacteria", "Gammaproteobacteria")) %>%
  group_by(Phylum_plot) %>%
  summarise(
    mean_abund = mean(Abundance),
    n_samples_present = sum(Abundance > 0),
    max_abund = max(Abundance)
  )

print("Alpha vs Gamma presence:")
print(alpha_gamma_check)
```


### Aggregate data per sample (sum over all Classes within each Phylum_plot category)
```python
plot_data <- class_data_plot %>%
  group_by(Sample, Phylum_plot, Clade) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Check dimensions
dim(plot_data)  
head(plot_data)
```

### Load phylogeny to get correct sample order
```python
rooted_red_sea_phylo <- read.tree("./Path/to/direcotory/Rooted_RAxML_RSHaplos_25.tre")

# Get phylogenetic order of samples
phylo_order <- rooted_red_sea_phylo$tip.label

# Factor Sample column in phylogenetic order
plot_data$Sample <- factor(plot_data$Sample, levels = phylo_order)

# Check for NAs (samples not in phylogeny)
sum(is.na(plot_data$Sample))
```

### Define colour palette
```python
phylum_colors <- c(
  "Alphaproteobacteria" = "#5E81AC",    # Blue grey
  "Gammaproteobacteria" = "#BF616A",    # Red
  "Cyanobacteria" = "#EBCB8B",          # Gold
  "Bacteroidota" = "#A3BE8C",           # Green
  "Chloroflexi" = "#B48EAD",            # Purple
  "Acidobacteriota" = "#88C0D0",        # Light blue
  "Actinobacteriota" = "#D08770",       # Orange/Red
  "Planctomycetota" = "#8FBCBB",        # Turquoise
  "Other" = "#4C566A"                   # Dark grey
)
```

### Scale abundance to 0-1 range
```python
plot_data_scaled <- plot_data %>%
  mutate(Abundance_scaled = Abundance / 100)  # Convert percentage to proportion

# Order taxa for legend - Custom order
taxa_custom_order <- c(
  "Alphaproteobacteria",    
  "Gammaproteobacteria",    
  "Cyanobacteria",
  "Bacteroidota",
  "Chloroflexi",
  "Acidobacteriota",
  "Actinobacteriota",
  "Planctomycetota",
  "Other"
)

# Factor with custom order
plot_data_scaled$Phylum_plot <- factor(plot_data_scaled$Phylum_plot, 
                                       levels = taxa_custom_order)

# Filter only major taxa
major_taxa_only <- plot_summary %>%
  filter(mean_abundance > 0.25, Phylum_plot != "Other") %>%
  pull(Phylum_plot)

print("Taxa to show in plot:")
print(major_taxa_only)

# Filter plot data. Keep major clades, remove NA clades
plot_data_filtered <- plot_data_scaled %>%
  filter(Phylum_plot %in% major_taxa_only) %>%
  filter(!is.na(Clade))  

# Check
sample_totals <- plot_data_filtered %>%
  group_by(Sample) %>%
  summarise(total = sum(Abundance_scaled))

print("Sample abundance range after filtering:")
print(paste("Min:", round(min(sample_totals$total), 3), "Max:", round(max(sample_totals$total), 3)))

# Update color palette (remove "Other")
phylum_colors_filtered <- phylum_colors[names(phylum_colors) %in% major_taxa_only]

# Order clades according to phylogeny
clade_phylo_order <- sample_info %>%
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
```python
Barplot_RS <- ggplot(plot_data_filtered, 
            aes(x = Sample, y = Abundance_scaled, fill = Phylum_plot)) +
  geom_bar(stat = "identity", 
           width = 1.0, 
           color = "white",
           linewidth = 0.1) +
  # FACET per clade - labels on top
  facet_grid(. ~ Clade, 
             scales = "free_x",      
             space = "free_x") +  
  # Colours
  scale_fill_manual(values = phylum_colors_filtered,
                    name = "Taxon") +
  # Y-axis
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  # Theme
  theme_minimal() +
  theme(
    # Axes
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    # Facet labels (clade names)
    strip.background = element_rect(fill = "grey", color = NA),
    strip.text = element_text(size = 6, face = "bold", color = "black", 
                              margin = margin(3, 0, 3, 0)),
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm"),
    # Panel background
    panel.spacing = unit(0.5, "lines"),  
    panel.grid = element_blank(),  # no grid
    panel.background = element_rect(fill = "grey92", color = NA),  
   # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    x = "Specimens (grouped by host clade)",
    y = "Relative abundance"
  )

# Display

print(Barplot_RS)
```








