rm(list = ls())
library("data.table")
library("dplyr")
library("metafor")
library("ggplot2")
library("ggrepel")

source('MicrobiomeAging_functions.R')

# Load all datasets into a named list
files <- list(
  FHSshot = "FHS/Results_species_FHS_shot.csv",
  MrOSshot = "MrOS/Results_species_MrOS_shot.csv",
  LLshot = "Lifelines_Frailty_results/Results_species_LL.csv",
  SOLshot = "SOL/Results_species_SOL_shot.csv"
)


dataframes <- load_and_prepare_data(
  files           = files,
  add_all = F
)

# Assign Mnumber
dataframes <- lapply(dataframes, mutate_summ_data)

# Load prefix lists
M3AgeDisplay <- data.table::fread("M3AgeDisplay.csv")$RISK
M3FIDisplay <- data.table::fread("M3FIDisplay.csv")$RISK

# Prepare valid prefixes
valid_age_prefixes <- gsub("^g__", "s__", M3AgeDisplay)
valid_age_prefixes <- gsub("genus", "species", valid_age_prefixes)

valid_fi_prefixes <- gsub("^g__", "s__", M3FIDisplay)
valid_fi_prefixes <- gsub("genus", "species", valid_fi_prefixes)

# Prepare risk levels for each outcome
risk_levels_age <- unique(unlist(lapply(dataframes, function(df) df$Variable)))
risk_levels_age <- risk_levels_age[sapply(risk_levels_age, function(risk) {
  any(startsWith(risk, valid_age_prefixes))
})]

risk_levels_fi <- unique(unlist(lapply(dataframes, function(df) df$Variable)))
risk_levels_fi <- risk_levels_fi[sapply(risk_levels_fi, function(risk) {
  any(startsWith(risk, valid_fi_prefixes))
})]

# Meta-analysis
meta_results_all <- run_meta_analysis(
  dataframes            = dataframes,
  outcome_levels        = c("age", "fi"),
  mnumber_levels        = c("m1", "m2", "m3", "m4"),
  population_levels     = "all",
  risk_levels_by_outcome = list(age = risk_levels_age, fi = risk_levels_fi),
  analysis_types        = "shot",
  ncohorts              = 4
) %>% subset(., !grepl('asv', RISK))


## Subset by outcome, population, and model
meta_results_all_m3 <- filter(meta_results_all, Mnumber == "m3")

# Load prefix lists from CSVs
M3AgeDisplay <- data.table::fread("M3AgeDisplay.csv")$RISK
M3FIDisplay <- data.table::fread("M3FIDisplay.csv")$RISK

# Process data
metaage <- process_data(meta_results_all, "age", "m3", M3AgeDisplay)
metafi <- process_data(meta_results_all, "fi", "m3", M3FIDisplay)

# List of outcomes to export
outcomes <- c("age", "fi", "cont", "mortality")

# List of valid prefixes
valid_prefixes <- list(
  age = M3AgeDisplay,
  fi = M3FIDisplay,
  cont = NULL,
  mortality = NULL
)




meta_results_age_m3 <- filter(metaage, Population == "all")
meta_results_fi_m3 <- filter(metafi, Population == "all")

#For validation
meta_results_age_val = filter(meta_results_age_m3, !is.na(delabel))
meta_results_fi_val  = filter(meta_results_fi_m3,  !is.na(delabel))
val_file = full_join(meta_results_age_val, meta_results_fi_val) %>%
  select(., c(RISK, OUTCOME2)) %>%
  rename(., 'Feature_Species'='RISK', 'Outcome'='OUTCOME2')
data.table::fwrite(val_file, "Updated_CHARGE_species.csv")


m3display <- full_join(meta_results_age_val, meta_results_fi_val) %>%
  build_m3display()

# Create separate displays for age and fi outcomes with both original and scaled CIs
m3agedisplayspecies = m3display %>% 
  filter(OUTCOME2 == 'age') %>%
  select(RISK, Studiesdir, nparticipants, B_CI, B_CI_original, pFDR_display, QE, QEp_display, I2, tau2, pFDR)

m3fidisplayspecies = m3display %>% 
  filter(OUTCOME2 == 'fi') %>%
  select(RISK, Studiesdir, nparticipants, B_CI, B_CI_original, pFDR_display, QE, QEp_display, I2, tau2, pFDR)

# Save updated files
write.csv(m3agedisplayspecies, "M3AgeDisplaySpecies.csv", row.names = FALSE)
write.csv(m3fidisplayspecies, "M3FIDisplaySpecies.csv", row.names = FALSE)

# List of outcomes to preview
outcomes <- c("age", "fi")

# Prefix lists per outcome
valid_prefixes <- list(
  age = M3AgeDisplay,
  fi = M3FIDisplay
)

# Save output
file_map <- list(
  age = "SuppTable6.csv",
  fi  = "SuppTable9.csv"
)

for (outcome in outcomes) {
  data_out <- meta_results_all %>%
    process_data_suppl(outcome = outcome, model = "m3") %>%
    filter(Population == "all") %>%
    prepare_supp_table() %>%
    mutate(across(where(is.list), unlist))
  
  if (nrow(data_out) > 0) {
    write_supp_table(data_out, file_map[[outcome]])
    cat("Saved:", file_map[[outcome]], "\n")
  } else {
    cat("No valid species found for outcome:", outcome, "\n")
  }
}

cohort_colours <- c(
  FHS = "#4477AA",         # blue
  MrOS = "#66CCEE",        # cyan
  LL = "#EE6677",          # red
  SOL = "#AA3377",         # purple
  "Meta-analysis" = "#BBBBBB" # grey
)

dataframes_all <- list(
  FHS  = dataframes$FHSshot  %>% filter(Datasplit == "all"),
  MrOS = dataframes$MrOSshot %>% filter(Datasplit == "all"),
  LL   = dataframes$LLshot   %>% filter(Datasplit == "all"),
  SOL  = dataframes$SOLshot  %>% filter(Datasplit == "all")
)

het_plots <- create_heterogeneity_forest_plots(
  meta_results    = full_join(meta_results_age_val, meta_results_fi_val) %>% subset(., pFDR < 0.05),
  dataframes_all  = dataframes_all,
  cohort_colours  = cohort_colours,
  save_patchwork  = FALSE,
  save_dir        = "../Paper/Supplements"
)
