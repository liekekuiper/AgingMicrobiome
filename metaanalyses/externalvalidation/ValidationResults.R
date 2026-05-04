rm(list = ls())
setwd('/Users/liekekuiper/Documents/Werk/VOILA/Microbiome/Meta analyse/SummaryStats/')

library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)
library(purrr)

source('../MicrobiomeAging_functions.R')

# -----------------------------------------------------------------------------
# VALIDATION-SPECIFIC HELPERS
# -----------------------------------------------------------------------------

#' Add pFDR, B_CI and formatted pFDR_display to a results data frame
add_stats <- function(df, p_col = "P", coef = "Coefficient") {
  df %>%
    mutate(
      pFDR        = p.adjust(.data[[p_col]], method = "fdr"),
      B_CI        = sprintf("%.2f (%.2f; %.2f)", .data[[coef]], LL, UL),
      pFDR_display = format_value(pFDR)   # from MicrobiomeAging_functions.R
    )
}

#' Filter for BMI model + significant M3 pairs
filter_m3 <- function(df, sig_pairs) {
  df %>%
    filter(str_detect(Model, "bmi")) %>%
    semi_join(sig_pairs, by = c("Variable" = "RISK", "Outcome"))
}

#' Full-join multiple datasets and label abundance threshold
combine_multi <- function(..., label) {
  reduce(list(...), full_join) %>%
    mutate(abprev = label) %>%
    add_stats()
}

#' Recode outcome strings to display labels
format_outcome <- function(x) {
  recode(x, age = "Chronological age", mortality = "All-cause Mortality")
}

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

genus_16    <- fread('FINRISK/R-version/16S/Results_FINRISK_16s_2025-11-14_g_subset.txt')
genus_me    <- fread('FINRISK/R-version/Metagenomics/Results_FINRISK_metagenomics_2025-09-26_g_subset.txt')
lenient_16s <- fread('FINRISK/R-version/Lenient/Results_FINRISK_16s_2025-11-04_g_subset.txt')
lenient_met <- fread('FINRISK/R-version/Lenient/Results_FINRISK_meta_2025-11-04_g_subset.txt')
simpson_16s <- fread('FINRISK/March2026/Results_FINRISK_new_additions_2026-03-17_simpson_mort_16s.txt')
simpson_met <- fread('FINRISK/March2026/Results_FINRISK_new_additions_2026-03-17_simpson_mort_Metagenome.txt')
species_me  <- fread('FINRISK/R-version/Metagenomics/Results_FINRISK_metagenomics_2025-09-26_s_subset.txt')

# -----------------------------------------------------------------------------
# SIGNIFICANT PAIRS FROM DISCOVERY
# -----------------------------------------------------------------------------

m3sig <- bind_rows(
  fread("M3AgeDisplay.csv") %>% mutate(Outcome = "age"),
  fread("M3FIDisplay.csv")  %>% mutate(Outcome = "mortality")
)
sig_pairs <- distinct(m3sig, RISK, Outcome)

m3spec <- bind_rows(
  fread("M3AgeDisplaySpecies.csv") %>% mutate(Outcome = "age",      QEp_display = as.character(QEp_display)),
  fread("M3FIDisplaySpecies.csv")  %>% mutate(Outcome = "mortality", QEp_display = as.character(QEp_display))
)
sig_pairs_spec <- distinct(m3spec, RISK, Outcome) %>%
  filter(!grepl("asv", RISK))

# -----------------------------------------------------------------------------
# GENUS VALIDATION
# pFDR corrected separately within strict and lenient (Option A)
# -----------------------------------------------------------------------------

g_strict <- combine_multi(
  filter_m3(genus_16,    sig_pairs),
  filter_m3(genus_me,    sig_pairs),
  filter_m3(simpson_16s, sig_pairs),
  filter_m3(simpson_met, sig_pairs),
  label = "strict"
)

missing_pairs <- anti_join(sig_pairs, g_strict, by = c("RISK" = "Variable", "Outcome"))

gl_add <- combine_multi(
  semi_join(filter_m3(lenient_16s, sig_pairs), missing_pairs, by = c("Variable" = "RISK", "Outcome")),
  semi_join(filter_m3(lenient_met, sig_pairs), missing_pairs, by = c("Variable" = "RISK", "Outcome")),
  label = "lenient"
)

g_both <- bind_rows(g_strict, gl_add)
fwrite(g_both, 'Validationgenus.csv')

# -----------------------------------------------------------------------------
# SPECIES VALIDATION
# -----------------------------------------------------------------------------

species_me_m3 <- species_me %>%
  filter_m3(sig_pairs_spec) %>%
  add_stats() %>%
  mutate(abprev = 'strict')

missing_spec <- anti_join(sig_pairs_spec, species_me_m3, by = c("RISK" = "Variable", "Outcome"))
if (nrow(missing_spec) > 0)
  message("Missing species combinations:\n",
          paste(missing_spec$RISK, missing_spec$Outcome, sep = " | ", collapse = "\n"))

fwrite(species_me_m3, 'ValidationSpecies.csv')

# -----------------------------------------------------------------------------
# BIOMARKERS
# -----------------------------------------------------------------------------

biomarkers <- bind_rows(
  transform(fread("FINRISK/March2026/FINRISK_16S_CoxPH_Comprehensive_withNcases.csv", data.table = FALSE), SeqMeth = "16S"),
  transform(fread("FINRISK/March2026/FINRISK_WGS_CoxPH_Comprehensive_withNcases.csv", data.table = FALSE), SeqMeth = "MG")
) %>%
  add_stats(p_col = "p.value", coef = "HR") %>%
  mutate(
    lrt_pfdr     = p.adjust(lrt_p, method = "fdr"),
    lrt_p_display = format_value(lrt_pfdr),
    Unit         = case_when(
      exposure == "genus_FI"    ~ "per percentage predicted frailty",
      exposure == "genus_FI_z"  ~ "per standard deviation",
      exposure == "AAgenus_FI"  ~ "per percentage predicted frailty (age-residual)",
      exposure == "AAgenus_FI_z"~ "per standard deviation (age-residual)"
    ),
    exposure2 = recode(exposure,
                       genus_FI     = "Raw genusFI",
                       genus_FI_z   = "Z-scored raw genusFI",
                       AAgenus_FI   = "Age-accelerated genusFI",
                       AAgenus_FI_z = "Z-scored age-accelerated genusFI"
    ),
    Model = gsub("raw_M|z_M|AA_M|AA_z_M", "Model ", model)
  )

# -----------------------------------------------------------------------------
# SUPPLEMENTARY TABLES
# -----------------------------------------------------------------------------

format_val_table <- function(df) {
  df %>%
    mutate(
      Variable2 = sapply(Variable, clean_label, table = TRUE),
      HR_CI     = if_else(Outcome == "mortality", B_CI, NA_character_),
      B_CI      = if_else(Outcome == "age",       B_CI, NA_character_),
      Outcome   = format_outcome(Outcome)
    )
}

# SuppTable5 — genus validation
g_both %>%
  arrange(Outcome, rev(abprev), pFDR) %>%
  format_val_table() %>%
  select(Outcome, abprev, Variable2, Variable, N, Ncases, B_CI, HR_CI, pFDR_display) %>%
  mutate(abprev = case_when(abprev == 'lenient' ~ 'less restrictive',
                            TRUE ~ abprev)) %>%
  rename(`Microbial Feature`        = Variable2,
         `Used abundance threshold` = abprev,
         `GreenGenes2 name`         = Variable,
         Nparticipants              = N,
         `Beta (CI)`                = B_CI,
         `Hazard Ratio (CI)`        = HR_CI,
         pFDR                       = pFDR_display) %>%
  write_supp_table("SuppTable5.csv")

# SuppTable7 — species validation
species_me_m3 %>%
  arrange(Outcome, pFDR) %>%
  format_val_table() %>%
  select(Outcome, Variable2, Variable, N, Ncases, B_CI, HR_CI, pFDR_display) %>%
  rename(`Microbial Feature` = Variable2,
         `GreenGenes2 name`  = Variable,
         Nparticipants       = N,
         `Beta (CI)`         = B_CI,
         `Hazard Ratio (CI)` = HR_CI,
         pFDR                = pFDR_display) %>%
  write_supp_table("SuppTable7.csv")

# SuppTable11 — biomarkers
biomarkers %>%
  select(SeqMeth, exposure2, Unit, Model, N, Ncases,
         B_CI, pFDR_display, concordance_full, lrt_p_display, c_index_diff) %>%
  rename(`Sequencing Method`                        = SeqMeth,
         Biomarker                                  = exposure2,
         `Hazard Ratio (CI)`                        = B_CI,
         Concordance                                = concordance_full,
         pFDR                                       = pFDR_display,
         `pFDR ANOVA Model without GenusFI`         = lrt_p_display,
         `ΔConcordance vs. model without GenusFI`   = c_index_diff) %>%
  write_supp_table("SuppTable11.csv")
