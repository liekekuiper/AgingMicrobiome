rm(list = ls())
library(data.table)
library(stringr)
library(dplyr)
library(openxlsx2)
library(ggplot2)
library(patchwork)
source('MicrobiomeAging_functions.R')

NOT_PRESENT <- "Not present at prevalence and abundance thresholds"
NCOLS     <- 14
DISC_COLS <- 1:8
MG_COLS   <- 9:11
S16_COLS  <- 12:14

# --------------------------------------------------------------------------
# LOAD & BUILD
# --------------------------------------------------------------------------
agegenus_CHARGE <- fread('M3AgeDisplay.csv')
agespeci_CHARGE <- fread('M3AgeDisplaySpecies.csv')
figenus_CHARGE  <- fread('M3FIDisplay.csv')
fispeci_CHARGE  <- fread('M3FIDisplaySpecies.csv')

val_genus     <- fread('Validationgenus.csv')
val_genus_16S <- subset(val_genus, N == 6606)
val_genus_MG  <- subset(val_genus, N == 6715)
val_species   <- fread('ValidationSpecies.csv')

age_gen  <- build_table(agegenus_CHARGE, val_genus_MG, val_genus_16S, "age", "age",      "Genus")
age_spec <- build_table(agespeci_CHARGE, val_species,  data.table(),  "age", "age",      "Species")
AgeRaw   <- apply_custom_sort(rbindlist(list(age_gen, age_spec), fill = TRUE))

fi_gen  <- build_table(figenus_CHARGE, val_genus_MG, val_genus_16S, "fi", "mortality", "Genus")
fi_spec <- build_table(fispeci_CHARGE, val_species,  data.table(),  "fi", "mortality", "Species")
FIRaw   <- apply_custom_sort(rbindlist(list(fi_gen, fi_spec), fill = TRUE))

rm(age_gen, age_spec, fi_gen, fi_spec)

# --------------------------------------------------------------------------
# PREPARE: asterisks, not-present flags, Nval column, column selection
# note: species rows have abprev_16S == NA (set in build_table)
#       genus "not present" rows have abprev_16S == NOT_PRESENT
# --------------------------------------------------------------------------
prep_for_output <- function(dt, is_fi = FALSE) {
  dt <- copy(dt)
  
  dt[, row_type := fifelse(is.na(abprev_16S), "species", "genus")]
  
  # Asterisk for lenient effect estimates
  dt[abprev_MG  == "lenient" & !is.na(B_CI_MG),  B_CI_MG  := paste0(B_CI_MG,  "*")]
  dt[abprev_16S == "lenient" & !is.na(B_CI_16S), B_CI_16S := paste0(B_CI_16S, "*")]
  
  # Not-present flags
  dt[, np_MG  := !is.na(abprev_MG)  & abprev_MG  == NOT_PRESENT]
  dt[, np_16S := !is.na(abprev_16S) & abprev_16S == NOT_PRESENT]
  
  # Build Nval display strings
  if (is_fi) {
    dt[, N_val_MG := fcase(
      np_MG  == TRUE,                       NOT_PRESENT,
      !is.na(N_MG) & !is.na(Ncases_MG),    paste0(N_MG,  " (", Ncases_MG,  ")"),
      default = NA_character_
    )]
    dt[, N_val_16S := fcase(
      np_16S == TRUE,                       NOT_PRESENT,
      !is.na(N_16S) & !is.na(Ncases_16S),  paste0(N_16S, " (", Ncases_16S, ")"),
      default = NA_character_
    )]
  } else {
    dt[, N_val_MG  := fcase(np_MG,  NOT_PRESENT, !is.na(N_MG),  as.character(N_MG),  default = NA_character_)]
    dt[, N_val_16S := fcase(np_16S, NOT_PRESENT, !is.na(N_16S), as.character(N_16S), default = NA_character_)]
  }
  
  # Clear effect/pFDR for not-present rows (will be merged in Excel)
  dt[np_MG  == TRUE, c("B_CI_MG",  "pFDR_display_MG")  := NA_character_]   
  dt[np_16S == TRUE, c("B_CI_16S", "pFDR_display_16S") := NA_character_]   
  
  dt[, .(
    row_type, np_MG, np_16S,
    c01 = sapply(Variable, clean_label, table = TRUE),   
    c02 = Studiesdir,
    c03 = nparticipants,
    c04 = B_CI,
    c05 = pFDR_display,
    c06 = QE,
    c07 = QEp_display,
    c08 = I2,
    c09 = N_val_MG,
    c10 = B_CI_MG,
    c11 = pFDR_display_MG,    
    c12 = N_val_16S,
    c13 = B_CI_16S,
    c14 = pFDR_display_16S   
  )]
}

age_prep <- prep_for_output(AgeRaw, is_fi = FALSE)
fi_prep  <- prep_for_output(FIRaw,  is_fi = TRUE)

# Insert Genus-level / Species-level section rows
add_sections <- function(dt) {
  blank <- function(label) {
    r           <- lapply(seq_len(ncol(dt)), function(i) NA_character_)
    names(r)    <- names(dt)
    r$row_type  <- "section"
    r$c01       <- label
    r$np_MG     <- FALSE
    r$np_16S    <- FALSE
    as.data.table(r)
  }
  rbindlist(list(
    blank("Genus-level"),
    dt[row_type == "genus"],
    blank("Species-level"),
    dt[row_type == "species"]
  ), fill = TRUE)
}

age_tbl <- add_sections(age_prep)
fi_tbl  <- add_sections(fi_prep)

# --------------------------------------------------------------------------
# RICH TEXT HEADER BUILDERS
# --------------------------------------------------------------------------

rt_N_sub <- function(sub) {
  fmt_txt("N", bold = FALSE) + fmt_txt(sub, bold = FALSE, vertAlign = "subscript")
}

rt_N_sub_parens <- function(sub1, sub2) {
  fmt_txt("N",   bold = FALSE) +
    fmt_txt(sub1, bold = FALSE, vertAlign = "subscript") +
    fmt_txt(" (N", bold = FALSE) +
    fmt_txt(sub2, bold = FALSE, vertAlign = "subscript") +
    fmt_txt(")",  bold = FALSE)
}

rt_pfdr <- function() {
  fmt_txt("p",   bold = FALSE, italic = TRUE) +
    fmt_txt("FDR", bold = FALSE)
}

rt_qep <- function() {
  fmt_txt("QE",  bold = FALSE) +
    fmt_txt("p",   bold = FALSE, italic = TRUE)
}

rt_i2 <- function() {
  fmt_txt("I",   bold = FALSE) +
    fmt_txt("2",   bold = FALSE, vertAlign = "superscript")
}

rt_plain <- function(x) fmt_txt(x, bold = FALSE)

# Column header lists for Age and FI
age_headers <- list(
  rt_plain("Microbial Feature"),                       #  1
  rt_N_sub_parens("studies", "dir"),                   #  2
  rt_N_sub("dis"),                                     #  3
  rt_plain("Beta (CI)"),                               #  4
  rt_pfdr(),                                           #  5
  rt_plain("QE"),                                      #  6
  rt_qep(),                                            #  7
  rt_i2(),                                             #  8
  rt_N_sub("val"),                                     #  9  MG
  rt_plain("Beta (CI)"),                               # 10
  rt_pfdr(),                                           # 11
  rt_N_sub("val"),                                     # 12  16S
  rt_plain("Beta (CI)"),                               # 13
  rt_pfdr()                                            # 14
)

fi_headers <- list(
  rt_plain("Microbial Feature"),
  rt_N_sub_parens("studies", "dir"),
  rt_N_sub("dis"),
  rt_plain("HR (CI)"),
  rt_pfdr(),
  rt_plain("QE"),
  rt_qep(),
  rt_i2(),
  rt_N_sub_parens("val", "cases"),                     #  9  MG
  rt_plain("HR (CI)"),
  rt_pfdr(),
  rt_N_sub_parens("val", "cases"),                     # 12  16S
  rt_plain("HR (CI)"),
  rt_pfdr()
)

# --------------------------------------------------------------------------
# RICH TEXT: VARIABLE LABELS (selective italic)
# --------------------------------------------------------------------------
rt_variable_label <- function(lbl) {
  if (is.na(lbl) || lbl == "") return(fmt_txt(""))
  
  # Diversity / uniqueness ~ plain
  if (grepl("^(Shannon|Simpson|Inverse Simpson|Chao1|Uniqueness)", lbl))
    return(fmt_txt(lbl))
  
  # "Genus sp. STRAIN ..."  italic genus word, roman remainder
  if (grepl("^\\S+ sp\\. ", lbl)) {
    genus <- regmatches(lbl, regexpr("^\\S+", lbl))
    rest  <- sub("^\\S+ ", "", lbl)
    return(fmt_txt(genus, italic = TRUE) + fmt_txt(paste0(" ", rest)))
  }
  
  # All other taxa ~ fully italic
  fmt_txt(lbl, italic = TRUE)
}

# --------------------------------------------------------------------------
# RICH TEXT: SCIENTIFIC NOTATION  "1.23x10^-04"  ~  subscript x, superscript exp
# Returns NULL when the value contains no scientific notation (write as-is).
# --------------------------------------------------------------------------
rt_sci <- function(v) {
  v <- as.character(v)
  if (!grepl("x10\\^", v)) return(NULL)
  
  # Strip leading zero from exponent: "^-04" ~ "^-4", "^03" ~ "^3"
  v <- gsub("\\^(-?)0+(\\d)", "^\\1\\2", v)
  
  # Split: prefix | x10^ | exponent | optional trailing "*"
  m <- regmatches(v, regexec("^(.*?)x10\\^(-?\\d+)(\\*?)$", v))[[1]]
  if (length(m) != 4) return(NULL)
  
  prefix   <- m[2]   # e.g. "1.23"
  exponent <- m[3]   # e.g. "-4"
  asterisk <- m[4]   # "*" or ""
  
  rt <- fmt_txt(prefix) +
    fmt_txt("x",      vertAlign = "subscript") +
    fmt_txt("10") +
    fmt_txt(exponent, vertAlign = "superscript")
  
  if (nchar(asterisk) > 0) rt <- rt + fmt_txt(asterisk)
  rt
}

# --------------------------------------------------------------------------
# EXCEL WRITER
# --------------------------------------------------------------------------
data_col_names <- paste0("c", sprintf("%02d", 1:NCOLS))

write_sheet <- function(wb, sheet_name, tbl, col_headers) {
  
  wb <- wb_add_worksheet(wb, sheet_name)
  
  section_rows <- which(tbl$row_type == "section")
  
  # --- ROW 1: group headers ---
  for (grp in list(list(DISC_COLS, "Discovery"),
                   list(MG_COLS,   "Metagenomics validation"),
                   list(S16_COLS,  "16S validation"))) {
    wb <- wb_merge_cells(wb, sheet_name, dims = wb_dims(rows = 1, cols = grp[[1]]))
    wb <- wb_add_data(wb, sheet_name, x = grp[[2]], dims = wb_dims(rows = 1, cols = min(grp[[1]])))
  }
  
  # --- ROW 2: column headers (rich text) ---
  for (i in seq_len(NCOLS))
    wb <- wb_add_data(wb, sheet_name, x = col_headers[[i]], dims = wb_dims(rows = 2, cols = i))
  
  # --- ROWS 3+: data ---
  for (ri in seq_len(nrow(tbl))) {
    sr  <- ri + 2
    row <- tbl[ri]
    
    if (row$row_type == "section") {
      wb <- wb_merge_cells(wb, sheet_name, dims = wb_dims(rows = sr, cols = 1:NCOLS))
      wb <- wb_add_data(wb, sheet_name, x = row$c01, dims = wb_dims(rows = sr, cols = 1))
      next
    }
    
    # CHANGED: column 1 — rich text with selective italic via rt_variable_label
    v1 <- row[["c01"]]
    if (!is.na(v1) && v1 != "NA")
      wb <- wb_add_data(wb, sheet_name,
                        x    = rt_variable_label(v1),
                        dims = wb_dims(rows = sr, cols = 1))
    
    # CHANGED: columns 2–14 — plain text or scientific-notation rich text
    for (ci in 2:NCOLS) {
      v <- row[[data_col_names[ci]]]
      if (!is.na(v) && v != "NA") {
        rt <- rt_sci(v)
        wb <- wb_add_data(wb, sheet_name,
                          x    = if (!is.null(rt)) rt else v,
                          dims = wb_dims(rows = sr, cols = ci))
      }
    }
    
    # Not-present: merge MG columns
    if (isTRUE(row$np_MG)) {
      wb <- wb_merge_cells(wb, sheet_name, dims = wb_dims(rows = sr, cols = MG_COLS))
      wb <- wb_add_data(wb, sheet_name, x = NOT_PRESENT, dims = wb_dims(rows = sr, cols = min(MG_COLS)))
    }
    
    # Not-present: merge 16S columns
    if (isTRUE(row$np_16S)) {
      wb <- wb_merge_cells(wb, sheet_name, dims = wb_dims(rows = sr, cols = S16_COLS))
      wb <- wb_add_data(wb, sheet_name, x = NOT_PRESENT, dims = wb_dims(rows = sr, cols = min(S16_COLS)))
      
      # CHANGED: species rows have no 16S data ~ show "-" in merged 16S block
    } else if (row$row_type == "species") {
      wb <- wb_merge_cells(wb, sheet_name, dims = wb_dims(rows = sr, cols = S16_COLS))
      wb <- wb_add_data(wb, sheet_name, x = "-", dims = wb_dims(rows = sr, cols = min(S16_COLS)))
    }
  }
  
  # --- STYLES ---
  total_rows <- nrow(tbl) + 2
  data_rows  <- setdiff(seq(3, total_rows), section_rows + 2)
  
  wb <- wb_add_border(wb, sheet_name,
                      dims       = wb_dims(rows = 1, cols = 1:NCOLS),
                      top_border = "thick", top_color = wb_color("black"))
  wb <- wb_add_font(wb, sheet_name,
                    dims = wb_dims(rows = 1, cols = 1:NCOLS), bold = "1")
  wb <- wb_add_cell_style(wb, sheet_name,
                          dims = wb_dims(rows = 1, cols = 1:NCOLS), horizontal = "center", vertical = "center")
  
  wb <- wb_add_border(wb, sheet_name,
                      dims          = wb_dims(rows = 2, cols = 1:NCOLS),
                      left_border   = NULL,
                      right_border  = NULL,
                      bottom_border = "thin", bottom_color = wb_color("black"))
  wb <- wb_add_cell_style(wb, sheet_name,
                          dims = wb_dims(rows = 2, cols = 1:NCOLS), horizontal = "center", vertical = "center")
  
  for (si in section_rows) {
    wb <- wb_add_font(wb, sheet_name,
                      dims = wb_dims(rows = si + 2, cols = 1:NCOLS), bold = "1")
    wb <- wb_add_border(wb, sheet_name,
                        dims       = wb_dims(rows = si + 2, cols = 1:NCOLS),
                        left_border   = NULL,
                        right_border  = NULL,
                        top_border = "thick", top_color = wb_color("black"))
  }
  
  
  wb <- wb_add_border(wb, sheet_name,
                      dims          = wb_dims(rows = total_rows, cols = 1:NCOLS),
                      top_border    = NULL, 
                      left_border   = NULL,
                      right_border  = NULL,
                      bottom_border = "thick", bottom_color = wb_color("black"))
  
  wb <- wb_set_col_widths(wb, sheet_name, cols = 1:NCOLS, widths = "auto")
  
  return(wb)
}

# --------------------------------------------------------------------------
# CREATE & SAVE — reassign wb at each step
# --------------------------------------------------------------------------
wb <- wb_workbook()
wb <- wb_set_base_font(wb, font_name = "Aptos Light", font_size = 10) 
wb <- write_sheet(wb, "Age Associations", age_tbl, age_headers)
wb <- write_sheet(wb, "FI Associations",  fi_tbl,  fi_headers)
wb_save(wb, "Microbiome_Associations_Validation.xlsx")
message("Formatted Tables 2 and 3 saved to Microbiome_Associations_Validation.xlsx")

### Heatmap 
# Heatmap data prep (bridges new Excel script and old heatmap code)
age_hm <- copy(AgeRaw)
fi_hm  <- copy(FIRaw)

# FI: combine N + Ncases into one column (mirrors old script)
fi_hm[, `N (Ncases) MG`  := fifelse(!is.na(N_MG)  & !is.na(Ncases_MG),
                                    paste0(N_MG,  " (", Ncases_MG,  ")"), NA_character_)]
fi_hm[, `N (Ncases) 16S` := fifelse(!is.na(N_16S) & !is.na(Ncases_16S),
                                    paste0(N_16S, " (", Ncases_16S, ")"), NA_character_)]
fi_hm[, c("N_MG", "Ncases_MG", "N_16S", "Ncases_16S") := NULL]

age_hm[] <- lapply(age_hm, as.character)
fi_hm[]  <- lapply(fi_hm,  as.character)

setnames(age_hm,
         old = c("Variable", "Studiesdir", "nparticipants", "B_CI", "pFDR_display",
                 "QE", "QEp_display", "I2",
                 "N_MG",  "B_CI_MG",  "pFDR_display_MG",  "abprev_MG",
                 "N_16S", "B_CI_16S", "pFDR_display_16S", "abprev_16S"),
         new = c("Microbial Feature", "Ncontributing studies (direction)",
                 "Nparticipants discovery", "Beta (95%-CI)",
                 "pFDR", "QE", "QEp", "I2",
                 "Nparticipants MG", "Beta (95%-CI) MG", "pFDR MG",
                 "Abundance threshold MG", "Nparticipants 16S",
                 "Beta (95%-CI) 16S", "pFDR 16S", "Abundance threshold 16S"))

setnames(fi_hm,
         old = c("Variable", "Studiesdir", "nparticipants", "B_CI", "pFDR_display",
                 "QE", "QEp_display", "I2",
                 "N (Ncases) MG",  "B_CI_MG",  "pFDR_display_MG",  "abprev_MG",
                 "N (Ncases) 16S", "B_CI_16S", "pFDR_display_16S", "abprev_16S"),
         new = c("Microbial Feature", "Ncontributing studies (direction)",
                 "Nparticipants discovery", "Beta (95%-CI)",
                 "pFDR", "QE", "QEp", "I2",
                 "N (Ncases) MG", "HR (95%-CI) MG", "pFDR MG",
                 "Abundance threshold MG", "N (Ncases) 16S",
                 "HR (95%-CI) 16S", "pFDR 16S", "Abundance threshold 16S"))

AgeAssociationsVal <- age_hm   
FIAssociationsVal  <- fi_hm
  
  
#  Rebuild data table 
age_dt <- as.data.table(AgeAssociationsVal)
fi_dt  <- as.data.table(FIAssociationsVal)
age_dt[, Outcome := "Age (years)"]
fi_dt[,  Outcome := "Frailty index (%)"]

dt <- rbindlist(list(age_dt, fi_dt), fill = TRUE)
dt[, pFDR_num := fifelse(pFDR == "<0.05", 0.045,
                         as.numeric(gsub("x10\\^", "e", pFDR)))]

#  Create Feature groups
dt[, FeatureGroup := dplyr::case_when(
  str_detect(`Microbial Feature`, "Shannon|Simpson|Chao1") ~ "\u03b1-diversity",
  str_detect(`Microbial Feature`, "Uniqueness")            ~ "Uniqueness",
  !grepl(" ", `Microbial Feature`) |
    grepl("^\\S+ sp\\. \\S+$", `Microbial Feature`)       ~ "Genus",
  TRUE                                                     ~ "Species"
)]

extract_effect <- function(x) as.numeric(str_extract(x, "^-?\\d+\\.?\\d*"))
dt[, Effect := extract_effect(`Beta (95%-CI)`)]
dt[FeatureGroup %in% c("\u03b1-diversity", 'Uniqueness'), Effect := Effect / 10]

# Add Significance & annotations
dt[, sig_MG         := is_sig(`pFDR MG`)]
dt[, sig_16S        := is_sig(`pFDR 16S`)]
dt[, any_sig        := sig_MG | sig_16S]
dt[, sig_MG_strict  := sig_MG  & (`Abundance threshold MG`  == "strict")]
dt[, sig_16S_strict := sig_16S & (`Abundance threshold 16S` == "strict")]

dt[, star   := sig_MG_strict & sig_16S_strict]
dt[FeatureGroup == "Species", star := sig_MG_strict]
# Shannon (Species) is MG-only — treat as species-level (star = MG strict only)
dt[`Microbial Feature` == "Shannon (Species)", star := sig_MG_strict]

dt[, dagger  := (sig_MG_strict) != (sig_16S_strict)]
dt[, ddagger := sig_16S & (`Abundance threshold 16S` == "lenient") &
     !sig_MG_strict & !sig_16S_strict]
dt[FeatureGroup == "Species",                              c("dagger","ddagger") := FALSE]
dt[`Microbial Feature` == "Shannon (Species)",             c("dagger","ddagger") := FALSE]

dt[, Annotation := ""]
dt[(star),                           Annotation := "*"]
dt[(!star) & (dagger),               Annotation := "\u2020"]
dt[(!star) & (!dagger) & (ddagger),  Annotation := "\u2021"]

# Add § to Shannon (Species) feature label to flag MG-only validation
dt[`Microbial Feature` == "Shannon (Species)", `Microbial Feature` := "Shannon (Species)\u00a7"]

# Sort 
dt[, GroupOrder := dplyr::case_when(
  FeatureGroup == "\u03b1-diversity" ~ 1,
  FeatureGroup == "Uniqueness"       ~ 2,
  FeatureGroup == "Genus"            ~ 3,
  FeatureGroup == "Species"          ~ 4,
  TRUE ~ 99L
)]

dt[, Nparticipants_discovery_num := as.numeric(`Nparticipants discovery`)]
dt <- dt[order(`Microbial Feature`, Outcome, -Nparticipants_discovery_num)]
dt <- dt[!duplicated(paste(`Microbial Feature`, Outcome))]
setorder(dt, GroupOrder, pFDR_num, any_sig, `Microbial Feature`)

plot_dt <- dt[, .(
  Feature = `Microbial Feature`,
  Outcome, Effect, Annotation, FeatureGroup, GroupOrder
)]


# Split by group & set factor levels 
# Relabel FeatureGroup
alpha_dt   <- set_feat_levels(plot_dt[FeatureGroup == "\u03b1-diversity"])
uniq_dt    <- set_feat_levels(plot_dt[FeatureGroup == "Uniqueness"])
genus_dt   <- set_feat_levels(plot_dt[FeatureGroup == "Genus"])
species_dt <- set_feat_levels(plot_dt[FeatureGroup == "Species"])

# Pad Uniqueness to two rows so a single feature doesn't get a giant tile
all_outcomes <- unique(dt$Outcome)
uniq_dt <- uniq_dt[, {
  have <- Outcome
  miss <- setdiff(all_outcomes, have)
  if (length(miss) == 0) .SD
  else rbind(.SD, .SD[1][, `:=`(Outcome = miss, Effect = NA_real_, Annotation = "")])
}, by = Feature]

# Colour scale 
fill_abs <- max(abs(plot_dt$Effect), na.rm = TRUE)

# Build panels 
p_alpha   <- make_panel(alpha_dt,   "#BBCCEE", fill_abs = fill_abs)
p_uniq    <- make_panel(uniq_dt,    "#CCEEFF", fill_abs = fill_abs)
p_genus   <- make_panel(genus_dt,   "#CCDDAA", show_x = TRUE, fill_abs = fill_abs)
p_species <- make_panel(species_dt, "#EEEEBB", show_x = TRUE, show_leg = TRUE, fill_abs = fill_abs)

# Assemble 
n_alpha   <- length(unique(alpha_dt$Feature))
n_uniq    <- length(unique(uniq_dt$Feature))
n_genus   <- length(unique(genus_dt$Feature))
n_species <- length(unique(species_dt$Feature))

left_rows  <- n_alpha + n_uniq + n_genus
right_rows <- n_species
pad_rows   <- max(0, left_rows - right_rows)

# Add a–d tags to each panel individually
p_alpha   <- p_alpha   + labs(tag = "a")
p_uniq    <- p_uniq    + labs(tag = "b")
p_genus   <- p_genus   + labs(tag = "c")
p_species <- p_species + labs(tag = "d")

# Style the tag in each panel's theme (add this line to make_panel, or patch here)
tag_theme <- theme(
  plot.tag          = element_text(size = 11, face = "bold", colour = "black"),
  plot.tag.position = "topleft"
)
p_alpha   <- p_alpha   + tag_theme
p_uniq    <- p_uniq    + tag_theme
p_genus   <- p_genus   + tag_theme
p_species <- p_species + tag_theme

p_left <- p_alpha / p_uniq / p_genus +
  plot_layout(heights = c(n_alpha, n_uniq, n_genus))

p_right <- p_species / plot_spacer() +
  plot_layout(heights = c(right_rows, pad_rows))

combined <- (p_left | p_right) +
  plot_layout(widths = c(1, 1))

#  Caption 
caption <- ggplot() +
  annotate(
    "text",
    x = 0, y = 1,
    hjust = 0, vjust = 1,
    size = 3, colour = "black",
    label = paste0(
      "* Genus-level: Validated in FINRISK in both MG and 16S datasets (strict threshold)\n",
      "* Species-level: Validated in FINRISK MG dataset (strict threshold);\n                          16S not considered for species-level features\n",
      "\u2020 Genus-level: Validated in one dataset (strict threshold)\n",
      "\u2021 Genus-level: 16S only, less restrictive threshold\n",
      "\u00a7 Species-level metric; validated in FINRISK MG only (16S not considered)"
    )
  ) +
  theme_void()

final_plot <- combined +
  inset_element(
    caption,
    left   = -2.1,
    bottom = 0.02,
    right  = 0.99,
    top    = 0.32,
    align_to = "panel"
  )

#  12. Save 
TILE_H     <- 0.48
total_rows <- max(left_rows, right_rows)
fig_height <- max(7, total_rows * TILE_H / 2.54 + 2.5)

ggsave("heatmap.png", final_plot,
       height = fig_height, width = 11, units = "in", dpi = 600, bg = "white")

message("Heatmap saved as heatmap.png")

