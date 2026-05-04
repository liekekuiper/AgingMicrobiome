rm(list=ls())
# Load required packages
library(data.table)
library(dplyr)
library(scales) 
library(ggplot2)
library(ggpattern)
library(sf)

# Define cohort name for files and for figure title
cohort_file <- c("MrOS", "FHS", "SOL")       # used in file names
cohort_label <- c("MrOS overlapping genera",
                  "FHS overlapping genera",
                  "HCHS/SOL overlapping genera") # used in figure title
basedir <- '/home/temp_mirna/l.m.kuiper/Microbiome/'

for(c in 1:length(cohort_file)){
  
  print(cohort_file[c])
  #setwd(paste0('~/Microbiome/', cohort_file, "/genus_export_", cohort_file, '/'))
  setwd(paste0(basedir, cohort_file[c], "/genus_export_", cohort_file[c], '/'))
  
  # Load and filter 16S and WGS prevalence data
  s16 <- fread(paste0("../", cohort_file[c], "_16S_Prev_Abundance_per_genus.csv"))
  s16 <- s16[`Percentage_with_≥0.01` >= 10, ]
  
  shot <- fread(paste0("../", cohort_file[c], "_WGS_Prev_Abundance_per_genus.csv"))
  shot <- shot[`Percentage_with_≥0.01` >= 10, ]
  
  # Load genus percentage table
  genus_perc <- fread(paste0(cohort_file[c], ".genus_table_percentages.tsv"))
  genus_perc$original_tax <- genus_perc$`#OTU ID`
  
  # Find common taxa between 16S and shotgun
  common_taxa <- intersect(s16$Original_Taxon, shot$Original_Taxon)
  
  # Keep only rows in genus_perc that are in the common taxa
  genus_perc_keep <- genus_perc[original_tax %in% common_taxa, ]
  
  cat('found at 16S:\n')
  print(sum(genus_perc$'16S' >0))
  
  cat('found at shotgun:\n')
  print(sum(genus_perc$'WGS' >0))
  
  
  # Print number shared genera
  cat("Shared genera:\n")
  print(nrow(genus_perc_keep))
  
  # Print number 16S genera
  cat("16S genera:\n")
  print(nrow(s16))
  
  # Print total percentage reads per platform
  cat("Percentage reads 16S:\n")
  print(sum(genus_perc_keep[[ "16S" ]]))
  
  # Print number shotgun genera
  cat("shotgun genera:\n")
  print(nrow(shot))
  
  cat("Percentage reads shotgun:\n")
  print(sum(genus_perc_keep[[ "WGS" ]]))
  
  # Create Venn diagram function
  create_venn <- function(data, file_name, main_title) {
    set1 <- data[[1]];  set2 <- data[[2]]
    n1   <- length(set1);  n2 <- length(set2)
    n12  <- length(intersect(set1, set2))
    n1_only <- n1 - n12;   n2_only <- n2 - n12
    
    col1 <- "#ffb000";  col2 <- "#648fff"
    
    # Build circle polygons with sf
    make_circle_sf <- function(cx, cy, r, n = 300) {
      theta <- seq(0, 2 * pi, length.out = n + 1)
      coords <- cbind(cx + r * cos(theta), cy + r * sin(theta))
      coords <- rbind(coords, coords[1, ])   # explicitly close the ring
      st_polygon(list(coords))
    }
    
    cx1 <- -0.38;  cx2 <- 0.38;  cy <- 0;  r <- 1
    c1 <- make_circle_sf(cx1, cy, r)
    c2 <- make_circle_sf(cx2, cy, r)
    
    only1  <- st_sfc(st_difference(c1, c2))
    only2  <- st_sfc(st_difference(c2, c1))
    inter  <- st_sfc(st_intersection(c1, c2))
    circ1  <- st_sfc(c1)
    circ2  <- st_sfc(c2)
    
    # Label positions
    lx1 <- cx1 - 0.8;  lx2 <- cx2 + 0.8
    
    p <- ggplot() +
      # Solid-filled exclusive regions
      geom_sf(data = only1, fill = col1, colour = NA) +
      geom_sf(data = only2, fill = col2, colour = NA) +
      # Overlap: col1 background + col2 diagonal stripes
      geom_sf_pattern(
        data            = inter,
        fill            = col1,
        colour          = NA,
        pattern         = "stripe",
        pattern_colour  = col2,
        pattern_fill    = col2,
        pattern_angle   = 45,
        pattern_density = 0.45,   # ~half stripe, half gap
        pattern_spacing = 0.06
      ) +
      # Circle outlines
      geom_sf(data = circ1, fill = NA, colour = col1, linewidth = 1.2) +
      geom_sf(data = circ2, fill = NA, colour = col2, linewidth = 1.2) +
      # Counts
      annotate("text", x = cx1 - 0.5, y = 0, label = n1_only, size = 12,
               fontface = "plain", family = "sans") +
      annotate("text", x = 0,          y = 0, label = n12,     size = 12,
               fontface = "plain", family = "sans") +
      annotate("text", x = cx2 + 0.5, y = 0, label = n2_only, size = 12,
               fontface = "plain", family = "sans") +
      ggtitle(main_title) +
      coord_sf(xlim = c(-1.8, 1.8), ylim = c(-1.2, 1.2), expand = FALSE) +
      theme_void() +
      theme(
        plot.title      = element_text(hjust = 0.5, size = 22,
                                       face = "plain", family = "sans"),
        plot.background = element_rect(fill = "white", colour = NA)
      )
    save(p, file = paste0('Venn_', cohort_file[c], '.rds'))
    ggsave(file_name, p, width = 5.5, height = 5.5, dpi = 1200)
  }
  
  # Generate Venn diagram
  create_venn(
    data = list("16S" = s16$Original_Taxon, "MG" = shot$Original_Taxon),
    file_name = paste0(cohort_file[c], "_venn_common_taxa_update.png"),
    main_title = cohort_label[c]
)
}
