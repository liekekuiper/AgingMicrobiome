rm(list=ls())
# Load required packages
install.packages('VennDiagram')
install.packages('scales') 
library(data.table)
library(dplyr)
library(VennDiagram)
library(scales) 

# Suppress VennDiagram logging
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Define cohort name for files and for figure title
cohort_file <- "FHS"       # used in file names
cohort_label <- "FHS"      # used in figure title

# Load and filter 16S and WGS prevalence data
s16 <- fread(paste0("../", cohort_file, "_16S_Prev_Abundance_per_genus.csv"))
s16 <- s16[`Percentage_with_≥0.01` >= 10, ]

shot <- fread(paste0("../", cohort_file, "_WGS_Prev_Abundance_per_genus.csv"))
shot <- shot[`Percentage_with_≥0.01` >= 10, ]

# Load genus percentage table
genus_perc <- fread(paste0(cohort_file, ".genus_table_percentages.tsv"))
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
  venn.diagram(
    x = data,
    imagetype = 'png',
    height = 5.5, width = 5.5, units = "in",
    category.names = c("16S", "MG"),
    col = c("#ffb000", "#648fff"),
    fill = c(alpha("#ffb000", 0.6), alpha("#648fff", 0.6)),
    cat.cex = 1, cat.fontfamily = "sans",
    cex = 4,
    main.fontfamily = "sans", fontfamily = "sans",
    main = main_title,
    main.cex = 5,
    filename = file_name,
    resolution = 1200
  )
}

# Generate Venn diagram
create_venn(
  data = list("16S" = s16$Original_Taxon, "MG" = shot$Original_Taxon),
  file_name = paste0(cohort_file, "_venn_common_taxa.png"),
  main_title = cohort_label
)
