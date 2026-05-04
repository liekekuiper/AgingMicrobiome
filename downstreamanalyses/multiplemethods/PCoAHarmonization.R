rm(list = ls())
library('qiime2R')
library('dplyr')
library('ggplot2')
library('ggtext')


# Define cohort name for files and for figure title

cohort_file <- c("MrOS", "FHS", "SOL")       # used in file names
cohort_label <- c("MrOS PCoA Weighted UniFrac",
                  "FHS PCoA Weighted UniFrac",
                  "HCHS/SOL PCoA Weighted UniFrac") # used in figure title
prep_qzv     <- 'weighted_prep_significance.qzv'
sample_qzv   <- 'weighted_sample_significance.qzv'
basedir <- '/home/temp_mirna/l.m.kuiper/Microbiome/'

for(c in 1:length(cohort_file)){
  print(cohort_file[c])
  setwd(paste0(basedir, cohort_file[c]))
  
  
  pcoa <- read_qza(paste0(cohort_file[c], '.asv.weighted.pc.qza'))
  
  dir.create("QZV_temp", showWarnings = FALSE)
  #system(paste("unzip -d QZV_temp/prep",   prep_qzv))
  #system(paste("unzip -d QZV_temp/sample", sample_qzv))
  
  adonis_prep   <- read.table(system("find ./QZV_temp/prep   -name adonis.tsv", intern = TRUE), header = TRUE, sep = "\t")
  adonis_sample <- read.table(system("find ./QZV_temp/sample -name adonis.tsv", intern = TRUE), header = TRUE, sep = "\t")
  
  # Convert PCoA results into a plottable data frame
  pcoa_df <- data.frame(pcoa$data$Vectors) %>%
    mutate(
      Preparation = case_when(
        grepl('16S', SampleID) ~ '16S',
        grepl('WGS', SampleID) ~ 'MG'
      )
    )
  
  # Pre-compute annotation positions
  x_range <- diff(range(pcoa_df$PC1))
  y_range <- diff(range(pcoa_df$PC2))
  
  ann_x  <- max(pcoa_df$PC1) - 0.015 * x_range  # 1.5% in from the right edge
  ann_y1 <- max(pcoa_df$PC2) - 0.05 * y_range   # 5% down from the top
  ann_y2 <- max(pcoa_df$PC2) - 0.20 * y_range   # 20% down from the top
  
  
  label_prep   <- paste0("<b>Sequencing method</b><br>R\u00B2: ", round(adonis_prep[1,   "R2"], 3))
  label_sample <- paste0("<b>Sample</b><br>R\u00B2: ", round(adonis_sample[1, "R2"], 3))
  
  plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Preparation, fill = Preparation)) +
    geom_point(
      shape = 21,
      size = 2,
      alpha = 0.8,
      stroke = 1
    ) +
    scale_color_manual(values = c('#ffb000', '#648fff')) +
    scale_fill_manual(values = c('#ffb000', '#648fff'), 
                      name = 'Sequencing method',
                      labels = c("16S rRNA gene amplicon", "shotgun metagenomics")) +
    labs(
      x     = paste0("PCo1 (", round(100 * pcoa$data$ProportionExplained[, 1], 2), "%)"),
      y     = paste0("PCo2 (", round(100 * pcoa$data$ProportionExplained[, 2], 2), "%)"),
      title = cohort_label  ) +
    theme_minimal(base_size = 12, base_family = "sans") +
    theme(
      plot.title      = element_text(face = "plain", hjust = 0.5, family = "sans"),
      legend.position = "bottom",
      legend.text     = element_text(size = 18, family = "sans"),  
      legend.title    = element_text(size = 18, family = "sans"),
      legend.key.size = unit(0.6, "cm")
    ) +
    guides(color = "none", fill = guide_legend(override.aes = list(stroke = 0, size = 10))) +
    annotate(
      "richtext",
      x           = ann_x,
      y           = ann_y1,
      label       = label_prep,
      size        = 5,
      lineheight  = 0.85,
      fill        = alpha("white", 0.55),
      label.color = NA,
      hjust       = 1
    ) +
    annotate(
      "richtext",
      x           = ann_x,
      y           = ann_y2,
      label       = label_sample,
      size        = 5,
      lineheight  = 0.85,
      fill        = alpha("white", 0.55),
      label.color = NA,
      hjust       = 1
    ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.3), axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
  
  save(plot, file = paste0(cohort_file[c], "_PCoA.rds"))
  ggsave(paste0(cohort_file[c], "_PCoA.png"), plot = plot,
         width = 6, height = 5, dpi = 900)
  
}
