rm(list=ls())
setwd("/Users/liekekuiper/Documents/Werk/VOILA/Microbiome/Meta analyse/SummaryStats/AGORA2/")
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(ggtext)
library(openxlsx)

# Load files
edt_full <- readxl::read_excel("ecological_direct_total_effects_M4_T0_1837.xlsx") %>%
  dplyr::select(
    VMHID,
    fullName,
    matches("Faecousia|vulgatus|Bacteroides_stercoris|Bacteroides_fragilis|Dysosmobacter_welbionis|Phocaeicola_plebeius|Blautia_wexlerae|Faecalibacterium_prausnitzii|Gemmiger_formicilis")
  )

#SCFA
scfa_ids <- c("ac", "but", "ppa", "isobut", "isoval")
scfa_data <- edt_full[edt_full$VMHID %in% scfa_ids, ]

#Secondary bile acids
sba_ids <- c("ddca", "HC02191", "12dhchol", "7dhcdchol")
sba_data <- edt_full[edt_full$VMHID %in% sba_ids, ]

#Other hypothesized toxic metabolites
other_ids <- c("indole", "tma", "h2s", "so3")
other_data <- edt_full[edt_full$VMHID %in% other_ids, ]

plot_forest <- function(data, dataset_name) {
  t_cols <- names(data)[grepl("^t_", names(data)) & !grepl("\\.lb$|\\.ub$", names(data))]
  species <- sub("^t_", "", t_cols)
  
  plot_df <- lapply(species, function(sp) {
    data %>%
      select(VMHID, fullName,
             Estimate = all_of(paste0("t_", sp)),
             LL = all_of(paste0("t_", sp, ".lb")),
             UL = all_of(paste0("t_", sp, ".ub"))) %>%
      mutate(Species = sp)
  }) %>% bind_rows()
  
  plot_df <- plot_df %>%
    mutate(
      fullName_clean = str_trim(str_extract(fullName, "^[^;]+")),
      Genus = str_extract(Species, "^[^_]+"),
      Species_part = str_extract(Species, "(?<=_).*")
    )
  
  plot_df <- plot_df %>%
    mutate(species_order = case_when(
      Species == "Faecalibacterium_prausnitzii" ~ 1,
      Species == "Gemmiger_formicilis" ~ 2,
      TRUE ~ 3
    )) %>%
    arrange(species_order, Species, fullName_clean) %>%
    mutate(y_axis = rev(seq_len(n())))
  
  # Assign genus/species labels to first and second row per species (from top to bottom)
  plot_df <- plot_df %>%
    group_by(Species) %>%
    arrange(desc(y_axis), .by_group = TRUE) %>%  # top to bottom
    mutate(
      Species_display = case_when(
        row_number() == 1 ~ paste0("*", Genus, "*"),
        row_number() == 2 ~ paste0("*", Species_part, "*"),
        TRUE ~ ""
      )
    ) %>%
    ungroup()
  
  # Horizontal line BELOW Gemmiger formicilis
  line_y <- plot_df %>%
    filter(Species == "Gemmiger_formicilis") %>%
    summarise(min_y = min(y_axis)) %>%
    pull(min_y) - 0.5
  
  ibm_colors <- c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000")
  
  ggplot(plot_df, aes(x = Estimate, y = y_axis, xmin = LL, xmax = UL, color = fullName_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = line_y, color = "black", linewidth = 0.6) +
    geom_point(size = 3.5) +
    geom_errorbarh(height = 0.2) +
    scale_y_continuous(
      breaks = plot_df$y_axis,
      labels = plot_df$Species_display,
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    scale_color_manual(values = ibm_colors) +
    labs(
      title = dataset_name,
      x = "Total effect on the community-level capacity\n(mmol/gram dry weight/hour)",
      y = NULL,
      color = "Metabolite"
    ) +
    guides(color = guide_legend(nrow = 2)) +
    theme_classic() +
    theme(
      axis.text.y = element_markdown(size = 15),
      axis.text.x = element_markdown(size = 12),
      axis.ticks.y = element_blank(), 
      axis.title.x = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, size = 18),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.box = "horizontal",
      legend.text = element_text(size = 15),      
      legend.title = element_text(size = 14, face = "bold")  
    )
}


# Example usage:
a = plot_forest(scfa_data, "Short chain fatty acids")
b = plot_forest(sba_data, "Secondary bile acids")
c = plot_forest(other_data, "Hypothesized toxic metabolites")

library('ggpubr')
#check = ggarrange(a, b, c, ncol = 3)
library(patchwork)
check <- (a | b | c) + plot_layout(ncol = 3)
ggsave('AGORA2.png',check, dpi = 900, units = 'in', width = 18, height = 8.5)

# Panels Fig4ab

# Extract display name from fullName column
extract_preferred_name <- function(fullname) {
  parts <- str_trim(unlist(str_split(fullname, "[;,]")))
  unwanted <- c("not degradable product", "undegradable product")
  
  # Check if first part is acceptable
  if (!any(str_detect(str_to_lower(parts[1]), unwanted))) {
    return(parts[1])
  }
  
  # Otherwise find first acceptable part
  for (part in parts[-1]) {
    if (!any(str_detect(str_to_lower(part), unwanted))) {
      return(part)
    }
  }
  
  # Fallback: return first part anyway
  return(parts[1])
}

# Vectorise it so it works with mutate()
extract_preferred_name_v <- Vectorize(extract_preferred_name)

plot_dumbbell <- function(df, bacteria_name, save_path = NULL) {
  
  t_col <- paste0("t_", bacteria_name)
  e_col <- paste0("e_", bacteria_name)
  d_col <- paste0("d_", bacteria_name)
  
  formatted_name <- str_replace_all(bacteria_name, "_", " ")
  
  plot_df <- df %>%
    filter(!is.na(.data[[t_col]]),
           !is.na(.data[[e_col]]),
           !is.na(.data[[d_col]])) %>%
    filter(abs(.data[[t_col]]) >= 1) %>%
    mutate(
      display_name = extract_preferred_name_v(fullName),
      direct       = .data[[d_col]],
      total        = .data[[t_col]],
      ecological   = .data[[e_col]]
    ) %>%
    arrange(total) %>%
    mutate(display_name = factor(display_name, levels = display_name))
  
  p <- ggplot(plot_df) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_segment(aes(x = direct, xend = total,
                     y = display_name, yend = display_name,
                     colour = "Ecological effect"),
                 linewidth = 1.2, lineend = "round") +
    geom_point(aes(x = direct, y = display_name, colour = "Direct effect"),
               size = 3.5) +
    geom_point(aes(x = total, y = display_name, colour = "Total effect"),
               size = 3.5) +
    scale_colour_manual(
      values = c(
        "Direct effect"     = "#ffb000",
        "Ecological effect" = "grey60",
        "Total effect"      = "#dc267f"
      ),
      guide = guide_legend(
        override.aes = list(
          shape     = c(16, NA, 16),
          linetype  = c(0,   1,  0),
          linewidth = c(0, 1.2,  0),
          size      = c(3.5, NA, 3.5)
        )
      )
    ) +
    labs(
      title  = paste0("*", formatted_name, "*"),
      x      = "Effect (mmol/gDW/hr)",
      y      = NULL,
      colour = NULL
    ) +
    theme_classic() +
    theme(
      plot.title         = element_markdown(hjust = 0.5, size = 18),
      axis.text.y        = element_text(size = 11),
      axis.text.x        = element_text(size = 12),
      axis.title.x       = element_text(size = 15),
      legend.position    = "bottom",
      legend.text        = element_text(size = 15),
      panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3)
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, dpi = 300, width = 8, height = 6, units = "in")
  }
  
  return(p)
}


# Species groupings
age_species <- c("Faecalibacterium_prausnitzii", "Gemmiger_formicilis")
fi_species   <- c("Phocaeicola_plebeius", "Phocaeicola_vulgatus",
                  "Bacteroides_fragilis", "Bacteroides_stercoris",
                  "Blautia_wexlerae")

select_species_cols <- function(df, species_vec) {
  pattern <- paste(species_vec, collapse = "|")
  df %>% select(VMHID, fullName, matches(pattern))
}

agora2_age <- select_species_cols(edt_full, age_species)
agora2_fi  <- select_species_cols(edt_full, fi_species)

p_fpr  <- plot_dumbbell(agora2_age, "Faecalibacterium_prausnitzii", 'Fprausnitzii_dumbbell.png')
p_gfo  <- plot_dumbbell(agora2_age, "Gemmiger_formicilis")
p_ppl  <- plot_dumbbell(agora2_fi,  "Phocaeicola_plebeius")
p_pvu  <- plot_dumbbell(agora2_fi,  "Phocaeicola_vulgatus")
p_bfr  <- plot_dumbbell(agora2_fi,  "Bacteroides_fragilis")
p_bst  <- plot_dumbbell(agora2_fi,  "Bacteroides_stercoris")
p_bwe  <- plot_dumbbell(agora2_fi,  "Blautia_wexlerae")

fig4 <- (p_fpr | p_gfo) / (a | b | c) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(heights = c(1, 1)) &
  theme(
    plot.tag = element_text(size = 20, face = "bold", colour = "black")
  )

ggsave(filename = 'fig4.png', fig4, units = 'in', height = 20, width = 18, dpi = 600)

# Excel sheet creation
edt_full <- edt_full %>%
  mutate(Category = case_when(
    VMHID %in% scfa_ids   ~ "Short chain fatty acids",
    VMHID %in% sba_ids    ~ "Secondary bile acids",
    VMHID %in% other_ids  ~ "Other hypothesized toxic metabolites",
    TRUE                  ~ "Other"
  ))

#Correlate ecological effect of age-associated species
eco_age = subset(edt_full, Category != 'Other')
cor(eco_age$e_Faecalibacterium_prausnitzii, 
    eco_age$e_Gemmiger_formicilis, method = 'spearman')

category_order <- c(
  "Short chain fatty acids",
  "Secondary bile acids",
  "Other hypothesized toxic metabolites",
  "Other"
)

supp_df <- edt_full %>%
  mutate(Metabolite = fullName) %>%
  select(-fullName) %>%
  pivot_longer(-c(VMHID, Metabolite, Category), names_to = "col", values_to = "value") %>%
  mutate(
    value   = as.numeric(value),
    effect  = str_extract(col, "^[edt](?=_)"),        # must be followed by _ 
    species = str_remove(col, "^[edt]_") %>%
      str_remove("\\.lb$") %>%
      str_remove("\\.ub$"),
    type = case_when(
      str_detect(col, "\\.lb$") ~ "lb",
      str_detect(col, "\\.ub$") ~ "ub",
      TRUE                      ~ "est"
    )
  ) %>%
  filter(!is.na(effect)) %>%                          # drop any rows that didn't match
  select(-col) %>%
  pivot_wider(
    names_from  = type,
    values_from = value,
    values_fn   = first                               # prevent list-cols from duplicates
  ) %>%
  mutate(
    cell_value = ifelse(
      is.na(est),
      "",
      paste0(round(est, 3), " (", round(lb, 3), "; ", round(ub, 3), ")")
    ),
    effect_label = case_when(
      effect == "d" ~ "Direct effect",
      effect == "e" ~ "Ecological effect",
      effect == "t" ~ "Total effect"
    ),
    col_name = paste0(str_replace_all(species, "_", " "), " ", effect_label)
  ) %>%
  select(VMHID, Category, Metabolite, col_name, cell_value) %>%
  pivot_wider(
    names_from  = col_name,
    values_from = cell_value,
    values_fn   = first
  ) %>%
  select(-VMHID) %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category, Metabolite) %>%
  mutate(Category = as.character(Category))



#Define order
all_species <- c(
  "Faecalibacterium_prausnitzii",
  "Gemmiger_formicilis",
  "Bacteroides_fragilis",
  "Bacteroides_stercoris",
  "Blautia_wexlerae",
  "Phocaeicola_plebeius",
  "Phocaeicola_vulgatus"
)

effect_labels <- c("Direct effect", "Ecological effect", "Total effect")

# Match the col_name format built during reshape: "Genus species Direct effect"
col_order <- c("Category", "Metabolite",
               as.vector(t(outer(str_replace_all(all_species, "_", " "), effect_labels,
                                 FUN = function(s, e) paste0(s, " ", e)))))

supp_df <- supp_df %>% select(any_of(col_order))

# Two header rows: row 1 = species spanners, row 2 = effect labels
header1 <- c("Category", "Metabolite",
             rep(str_replace_all(all_species, "_", " "), each = 3))
header2 <- c("", "",
             rep(effect_labels, length(all_species)))

wb <- createWorkbook()
addWorksheet(wb, "Supplement")

style_header_species <- createStyle(
  fontName       = "Aptos Light", fontSize = 10, fontColour = "#FFFFFF",
  fgFill         = "#2C3E50", halign = "CENTER", valign = "CENTER",
  textDecoration = "bold", border = "Bottom", borderColour = "#FFFFFF",
  wrapText       = TRUE
)
style_header_effect <- createStyle(
  fontName       = "Aptos Light", fontSize = 10, fontColour = "#FFFFFF",
  fgFill         = "#566573", halign = "CENTER", valign = "CENTER",
  textDecoration = "bold", border = "Bottom", borderColour = "#FFFFFF"
)

style_header_species <- createStyle(
  fontName       = "Aptos Light", fontSize = 10, fontColour = "#000000",
  fgFill         = "#FFFFFF", halign = "CENTER", valign = "CENTER",
  textDecoration = "bold", border = "Bottom", borderColour = "#000000",
  wrapText       = TRUE
)
style_header_effect <- createStyle(
  fontName       = "Aptos Light", fontSize = 10, fontColour = "#000000",
  fgFill         = "#FFFFFF", halign = "CENTER", valign = "CENTER",
  textDecoration = "bold", border = "Bottom", borderColour = "#000000"
)

style_metabolite <- createStyle(
  fontName       = "Aptos Light", fontSize = 10,
  textDecoration = "italic", halign = "LEFT"
)
style_data <- createStyle(
  fontName = "Aptos Light", fontSize = 10, halign = "CENTER"
)
style_stripe <- createStyle(
  fontName = "Aptos Light", fontSize = 10,
  halign   = "CENTER", fgFill = "#F2F3F4"
)

writeData(wb, "Supplement", as.data.frame(t(header1)), startRow = 1, colNames = FALSE)
writeData(wb, "Supplement", as.data.frame(t(header2)), startRow = 2, colNames = FALSE)
writeData(wb, "Supplement", supp_df,                   startRow = 3, colNames = FALSE)

n_rows    <- nrow(supp_df)
n_cols    <- ncol(supp_df)
data_rows <- 3:(3 + n_rows - 1)

addStyle(wb, "Supplement", style_header_species, rows = 1,        cols = 1:n_cols,  gridExpand = TRUE)
addStyle(wb, "Supplement", style_header_effect,  rows = 2,        cols = 1:n_cols,  gridExpand = TRUE)
addStyle(wb, "Supplement", style_metabolite,     rows = data_rows, cols = 1:2,       gridExpand = TRUE)
addStyle(wb, "Supplement", style_data,           rows = data_rows, cols = 3:n_cols,  gridExpand = TRUE)

#stripe_rows <- data_rows[seq(2, length(data_rows), 2)]
#addStyle(wb, "Supplement", style_stripe, rows = stripe_rows, cols = 3:n_cols, gridExpand = TRUE)

# After your existing addStyle calls, add this:
style_border_top <- createStyle(
  fontName       = "Aptos Light", fontSize = 10, fontColour = "#000000",
  fgFill         = "#FFFFFF", halign = "CENTER", valign = "CENTER",
  textDecoration = "bold",
  border         = "Top", borderColour = "#000000"   # ← top border on row 2
)
addStyle(wb, "Supplement", style_border_top, rows = 2, cols = 1:n_cols, 
         gridExpand = TRUE, stack = TRUE)             # stack = TRUE preserves existing styles

# Spanners start at col 3 (Category=1, Metabolite=2)
for (i in seq_along(all_species)) {
  start_col <- 3 + (i - 1) * 3
  mergeCells(wb, "Supplement", cols = start_col:(start_col + 2), rows = 1)
}

setColWidths(wb, "Supplement", cols = 1,        widths = 32)
setColWidths(wb, "Supplement", cols = 2,        widths = 40)
setColWidths(wb, "Supplement", cols = 3:n_cols, widths = 22)

freezePane(wb, "Supplement", firstActiveRow = 3)
saveWorkbook(wb, "../../Paper/Supplements/SuppTable10.xlsx", overwrite = TRUE)

