rm(list=ls())

library(ggplot2)
library(dplyr)
library(readr)

setwd('/home/temp_mirna/l.m.kuiper/Microbiome')

# Load the scores from MrOS, SOL, Framingham Heart Study, and Rotterdam Study
scores_server1 <- read_csv("combined_pca_scores.csv") %>%
  mutate(Cohort = case_when(Cohort == 'RS' ~ 'RS 16S',
                            TRUE ~ Cohort))

# Load the scores from the Doetinchem Cohort Study
scores_dcs <- read_csv("pca_scores_DCS_common_genera.csv") %>%
  mutate(Cohort = 'DCS 16S')

# Load the scores from LifeLines
scores_ll <- read_csv("pca_scores_DAGIII_common_genera.csv") %>%
  mutate(Cohort = "LL MG")

# Combine scores into one dataframe (outer joins)
scores_combined <- scores_server1 %>%
  full_join(scores_dcs) %>%
  full_join(scores_ll)

# Reshuffle rows to avoid order bias
set.seed(123)  # optional for reproducibility
scores <- scores_combined %>%
  slice_sample(prop = 1)

# Load explained variance per PCA
explained_variance <- read_csv("explained_variance.csv")

pc1_var <- explained_variance$`0`[1] * 100
pc2_var <- explained_variance$`0`[2] * 100

# Define cohort colors (matching your upset_query palette)
cohort_colors <- c(
  "FHS 16S" = "#cc6677",
  "FHS MG"  = "#332288",
  "MrOS 16S"= "#ddcc77",
  "MrOS MG" = "#117733",
  "RS 16S"  = "#88ccee",
  "LL MG"   = "#882255",
  "DCS 16S" = "#44AA99",
  "HCHS/SOL 16S" = "#999933",
  "HCHS/SOL MG"  = "#AA4499"
)

cohort_colors <- c(
  "FHS 16S" = "#882255",
  "FHS MG"  = "#AA4499",
  "MrOS 16S"= "#332288",
  "MrOS MG" = "#88ccee",
  "RS 16S"  = "#117733",
  "LL MG"   = "#cc6677",
  "DCS 16S" = "#44AA99",
  "HCHS/SOL 16S" = "#999933",
  "HCHS/SOL MG"  = "#ddcc77"
)

# Ensure alphabetical legend order
scores$Cohort <- factor(scores$Cohort,
                        levels = sort(unique(scores$Cohort)))

# Plot PCA
p <- ggplot(scores, aes(x = PC1, y = PC2, color = Cohort, fill = Cohort)) +
  geom_point(
    shape = 21,
    size = 2,
    alpha = 0.8,
    stroke = 1
  ) +
  scale_color_manual(values = cohort_colors, guide = "none") +
  scale_fill_manual(values = cohort_colors) +
  labs(
    title = "PCA of CLR-Transformed Data for Top 7 Genera (Grouped by Cohort)",
    x = sprintf("PC1 (%.1f%%)", pc1_var),
    y = sprintf("PC2 (%.1f%%)", pc2_var),
    fill = NULL   # removes legend title
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      override.aes = list(
        color = NA,   # removes outline in legend
        stroke = 0,    # ensures no border
        size = 6
      )
    )
  ) +
  theme_minimal(base_size = 22, base_family = "sans") +
  theme(
    plot.title = element_text(
      face = "plain",
      hjust = 0.5,
      family = "sans"
    ),
    legend.position = c(0.85, 0.9),
    #legend.background = element_rect(
    #  fill = "white",
    #  color = "black",
    #  linewidth = 0.3
    #),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Save plot
ggsave(
  filename = "combined_pca_plot.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 600
)

print(p)
save(p, file = 'combined_pca_plot.rds')
