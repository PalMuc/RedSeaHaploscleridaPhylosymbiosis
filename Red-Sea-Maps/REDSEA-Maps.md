# R-Script for configuring Red Sea Maps

### Setup libraries
```python
library(tidyverse)
library(readxl)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(cowplot)
```

### Load the data
```python
sample_data <- read_excel("C:/Path/to/directory/R-Phylosymbiosis/Supplementary_Table_1.xlsx")

# Fix coordinate issues
sample_data <- sample_data %>%
  mutate(Longitude = ifelse(Longitude > 100, Longitude / 1000000, Longitude))

# Remove any samples with missing or invalid coordinates
sample_data_clean <- sample_data %>%
  filter(!is.na(Latitude), !is.na(Longitude),
         Longitude >= 30, Longitude <= 50,
         Latitude >= 10, Latitude <= 32)
```

### Create subsets
```python
# Map A: G15 + G16 + G17
map_A_data <- sample_data_clean %>%
  filter(Clade %in% c("G15", "G16", "G17"))

# Map B: G04 + G07 + G19 
map_B_data <- sample_data_clean %>%
  filter(Clade %in% c("G04","G07","G19"))

# Map C: G01 + G06 + G12 
map_C_data <- sample_data_clean %>%
  filter(Clade %in% c("G01","G06","G12"))
```

### Prepare spatial data
```python
# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Convert to sf objects
map_A_sf <- st_as_sf(map_A_data, 
                     coords = c("Longitude", "Latitude"), 
                     crs = 4326, agr = "constant")

map_B_sf <- st_as_sf(map_B_data, 
                     coords = c("Longitude", "Latitude"), 
                     crs = 4326, agr = "constant")

map_C_sf <- st_as_sf(map_C_data, 
                     coords = c("Longitude", "Latitude"), 
                     crs = 4326, agr = "constant")
                     
# Define Red Sea bounding box
xlim_red_sea <- c(32, 44)
ylim_red_sea <- c(12, 30)
```

### Define colours
```python
# Colors for each clade
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

# Alternative: Single color per map (cleaner look)
map_A_colour <- c("G15" = "midnightblue", "G16" = "lightcyan3", "G17" = "firebrick2")  
map_B_colour <- c("G04" = "purple3", "G07" = "cornflowerblue", "G19" = "darkseagreen")  
map_C_colour <- c("G01" = "blue3", "G06" = "greenyellow", "G12" = "cyan")
```

## Map A

```python
map_A <- ggplot(data = world) +
  geom_sf(fill = "grey90", color = "grey60", linewidth = 0.3) +
  geom_sf(data = map_A_sf, 
          aes(fill = Clade),
          size = 3.5, 
          shape = 21, 
          color = "black",
          stroke = 0.8,
          alpha = 0.8) +
  scale_fill_manual(values = map_A_colour,
                    name = "Host clade",
                    labels = c("G15", 
                               "G16",
                               "G17")) +
  coord_sf(xlim = xlim_red_sea, ylim = ylim_red_sea, expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.3,
                   text_cex = 0.8, line_width = 1) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  labs(title = "A. Clades G15 + G16 + G17",
       subtitle = paste0("n = ", 
                         nrow(map_A_data))) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
    panel.background = element_rect(fill = "aliceblue"),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

print(map_A)
```

## Map B
```python
map_B <- ggplot(data = world) +
  geom_sf(fill = "grey90", color = "grey60", linewidth = 0.3) +
  geom_sf(data = map_B_sf, 
          aes(fill = Clade),
          size = 3.5, 
          shape = 21, 
          color = "black",
          stroke = 0.8,
          alpha = 0.8) +
  scale_fill_manual(values = map_B_colour,
                    name = "Host clade",
                    labels = c("G04", 
                               "G07",
                               "G19")) +
  coord_sf(xlim = xlim_red_sea, ylim = ylim_red_sea, expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.3,
                   text_cex = 0.8, line_width = 1) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  labs(title = "B. Clades G04 + G07 + G19",
       subtitle = paste0("n = ", 
                         nrow(map_B_data))) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
    panel.background = element_rect(fill = "aliceblue"),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

print(map_B)
```

## Map C
```python
map_C <- ggplot(data = world) +
  geom_sf(fill = "grey90", color = "grey60", linewidth = 0.3) +
  geom_sf(data = map_C_sf, 
          aes(fill = Clade),
          size = 3.5, 
          shape = 21, 
          color = "black",
          stroke = 0.8,
          alpha = 0.8) +
  scale_fill_manual(values = map_C_colour,
                    name = "Host clade",
                    labels = c("G01", 
                               "G06",
                               "G12")) +
  coord_sf(xlim = xlim_red_sea, ylim = ylim_red_sea, expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.3,
                   text_cex = 0.8, line_width = 1) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  labs(title = "B. Clades G01 + G06 + G12",
       subtitle = paste0("n = ", 
                         nrow(map_C_data))) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
    panel.background = element_rect(fill = "aliceblue"),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

print(map_C)
```

### For combined figure:
```python
combined_horizontal <- plot_grid(
  map_A + theme(legend.position = "none"),
  map_B + theme(legend.position = "none"),
  map_C + theme(legend.position = "none"),
  ncol = 3,
  labels = c("", "", ""),
  align = "h",
  axis = "tb"
)

# Extract legends
legend_A <- get_legend(map_A)
legend_B <- get_legend(map_B)
legend_C <- get_legend(map_C)

# Combine legends
combined_legends <- plot_grid(legend_A, legend_B, legend_C, ncol = 3)

# Final combined plot with legends below
combined_with_legends <- plot_grid(
  combined_horizontal,
  combined_legends,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

ggsave("RedSea_Maps_Combined.pdf", combined_with_legends, 
       width = 18, height = 8)
ggsave("RedSea_Maps_Combined.svg", combined_with_legends, 
       width = 18, height = 8, dpi = 300)
```


