rm(list = ls())
library("data.table")
library("dplyr")
library("ggrepel")
library("tidyr")
library("ComplexUpset")
library('stringr')
library('metafor')
library('patchwork')

source('MicrobiomeAging_functions.R')

files <- list(
  FHS16s   = "FHS/Results_FHS_16s.csv",
  FHSshot  = "FHS/Results_FHS_shot.csv",
  MrOS16s  = "MrOS/Results_MrOS_16s.csv",
  MrOSshot = "MrOS/Results_MrOS_shot.csv",
  RS16s    = "RS/Results_RS_16s.csv",
  DCS16s   = "DCS/Results_DCS_16s.csv",
  LLshot   = "Lifelines_Frailty_results/intermediatefiles/Results_age_5_DAG3_2026-02-13.csv",
  SOL16s   = "SOL/Results_SOL_16s.csv",
  SOLshot  = "SOL/Results_SOL_shot.csv"
)

mean_cohort_age = list(
  FHS16s = 55.6,
  FHSshot = 55.7,
  MrOS16s = 84.2,
  MrOSshot = 84.1,
  RS16s = 62.6,
  DCS16s = 67.2,
  LLshot = 51.5,
  SOL16s = 57.1,
  SOLshot = 56.7
)

dataframes <- load_and_prepare_data(
  files           = files,
  expand_fn       = expand_MrOS,
  expand_datasets = c("MrOS16s", "MrOSshot"),
  mean_cohort_age_list = mean_cohort_age
)

dataframes <- lapply(dataframes, mutate_summ_data)

##### Random effect meta analyses ####
re_meta_results <- run_meta_analysis(
  dataframes       = dataframes,
  outcome_levels   = c("age", "fi"),
  mnumber_levels   = c("m1", "m2", "m3", "m4"),
  population_levels = unique(dataframes$LLshot$Datasplit),
  analysis_types   = c("16s", "shot", "all")
)

meta_results_df = subset(re_meta_results, Population == 'all' & !grepl('species', RISK))
# FDR adjust
meta_results <- meta_results_df %>%
  mutate(pFDR = p.adjust(P, method = "fdr"))

fwrite(meta_results, "meta_results_genus.csv")
m3_meta_results_df = subset(meta_results, Mnumber == 'm3' & Analysis_Type == "all" & Population == 'all')
m3sig = subset(meta_results, Mnumber == 'm3' & Analysis_Type == "all" & Population == 'all' & pFDR < 0.05 & ndirection > 1)


# sig combinations m3
sig_pairs <- unique(m3sig[c("RISK", "OUTCOME2")])

# filter m4
m4sig <- subset(
  meta_results,
  Mnumber == "m4" &
    Analysis_Type == "all" &
    Population == "all" &
    pFDR < 0.05 &
    interaction(RISK, OUTCOME2) %in% interaction(sig_pairs$RISK, sig_pairs$OUTCOME2)
)

supptabl2 <- m3_meta_results_df %>% process_data_suppl(outcome = "age") %>% arrange(pFDR) %>% prepare_supp_table()
supptabl8 <- m3_meta_results_df %>% process_data_suppl(outcome = "fi")  %>% arrange(pFDR) %>% prepare_supp_table()
dietsens  <- meta_results %>%
  filter(Mnumber == 'm4', Analysis_Type == 'all',
         paste(RISK, OUTCOME2) %in% paste(m3sig$RISK, m3sig$OUTCOME2)) %>%
  process_data_suppl() %>%
  prepare_supp_table()

write_supp_table(supptabl2, "SuppTable2.csv")
write_supp_table(supptabl8, "SuppTable8.csv")
write_supp_table(dietsens,  "SuppTable3.csv")

#Create Display version for Model 3 - Main model
m3display <- m3sig %>% build_m3display()

# Save feature-genus mapping
val_genus = select(m3display, c(RISK, OUTCOME2)) %>%
  rename(., 'Feature_Genus'='RISK', 'Outcome'='OUTCOME2')
data.table::fwrite(val_genus, "Updated_CHARGE_genus.csv")

# Create separate displays for age and fi outcomes with both original and scaled CIs
m3agedisplay = m3display %>% 
  filter(OUTCOME2 == 'age') %>%
  select(Variable,RISK, Studiesdir, nparticipants, B_CI, B_CI_original, pFDR_display, QE, QEp_display, I2, tau2, pFDR)

m3betafidisplay = m3display %>% 
  filter(OUTCOME2 == 'fi') %>%
  select(Variable,RISK, Studiesdir, nparticipants, B_original, B_CI, B_CI_original, pFDR_display, QE, QEp_display, I2, tau2, pFDR) %>%
  rename(., B = B_original)
m3fidisplay = select(m3betafidisplay, -B)

# Save updated files
write.csv(m3agedisplay, "M3AgeDisplay.csv", row.names = FALSE)
write.csv(m3fidisplay, "M3FIDisplay.csv", row.names = FALSE)
write.csv(m3betafidisplay, "M3_B_FIDisplay.csv", row.names = FALSE)
#### Subgroup analyses ####
#### Heterogeneity

cohort_names <- c("FHS", "MrOS", "RS", "DCS", "LL", "SOL", "Meta-analysis")

cohort_colours <- c(
  FHS = "#4477AA",         # blue
  MrOS = "#66CCEE",        # cyan
  RS = "#228833",          # green
  DCS = "#CCBB44",         # yellow
  LL = "#EE6677",          # red
  SOL = "#AA3377",         # purple
  "Meta-analysis" = "#BBBBBB" # grey
)

# Define cohorts used in meta-analyses
dataframes_all <- list(FHS = dataframes$FHSall, 
                       MrOS = dataframes$MrOSall, 
                       RS = dataframes$RSall, 
                       DCS = dataframes$DCSall, 
                       LL = dataframes$LLall, 
                       SOL = dataframes$SOLall)

het_plots <- create_heterogeneity_forest_plots(
  meta_results         = m3sig,
  dataframes_all       = dataframes_all,
  cohort_colours       = cohort_colours,
  save_dir             = "../Paper/Supplements/"
)


meta_fdr_clean <- run_meta_regression(
  sig_pairs      = m3sig[, c("RISK", "OUTCOME2")],
  dataframes_all = dataframes_all
)

supptable4 = meta_fdr_clean %>%
  mutate(
    Variable = sapply(RISK, clean_label, TRUE), 
    Comparison = case_when(Population_Group == 'women' ~ 'Women (reference)',
                           Population_Group == 'men' ~ "Men, reference: women",
                           Population_Group == 'age_3' ~ 'Age 50-60 (reference)',
                           Population_Group == 'age_1' ~ 'Age < 40 years, reference: Age 50-60',
                           Population_Group == 'age_2' ~ 'Age 40-50, reference: Age 50-60',
                           Population_Group == 'age_4' ~ 'Age 60-70, reference: Age 50-60',
                           Population_Group == 'age_5' ~ 'Age 70+, reference: Age 50-60'),
    # Confidence intervals for original values
    B_CI = case_when(OUTCOME2 == 'fi' ~ paste0(sprintf(B, fmt = "%.4f"), " (", sprintf(LL, fmt = "%.4f"), "; ", sprintf(UL, fmt = "%.4f"), ")"),
                     TRUE ~ paste0(sprintf(B, fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), "; ", sprintf(UL, fmt = "%.2f"), ")")),
    # Format p-values
    pFDR_display = format_value(pFDR),
    QM_display = case_when(is.na(QM_display) ~ NA_character_,
                           TRUE ~ sprintf(QM_display, fmt = "%.2f")),
    QE_display = case_when(is.na(QE_display) ~ NA_character_, 
                           TRUE ~ sprintf(QE_display, fmt = "%.2f")),
    QEp_display = format_value(QEpFDR_display),
    QMpFDR_display = format_value(QMpFDR_display),
    I2 = sprintf(I2, fmt= "%.2f"),
    tau2 = sprintf(tau2, fmt= "%.2f"),
    Outcome = case_when(OUTCOME2 == 'age' ~ 'Chronological age',
                        OUTCOME2 == 'fi' ~ 'Frailty index'),
    Studiesdir = paste0(nstudies, " (", ndirection, ")")
  )  %>%
  select(Outcome, Comparison, Variable, RISK, Studiesdir, nparticipants, B_CI, pFDR_display, QM_display, QMpFDR_display, R2_display, QE_display, QEp_display, I2, tau2) %>%
  rename(., 
         `Microbial Feature` = Variable,
         `GreenGenes2 name` = RISK,
         `Nstudies (Ndirection)`= Studiesdir,
         Nparticipants = nparticipants,
         `Beta (CI)` = B_CI,
         pFDR = pFDR_display,
         QM = QM_display,
        `QM pFDR` = QMpFDR_display,
        R2 = R2_display,
        QE = QE_display,
         `QE p-value` =QEp_display) 
write_supp_table(supptable4, "SuppTable4.csv")

#### Figures ####
## Make figure sequencing type
# Generate and save plots for 'age' and 'fi'
forest_plot_age <- create_forest_plot(meta_results, m3sig, "age", "Chronological age")
forest_plot_fi <- create_forest_plot(meta_results, m3sig, "fi", "the Frailty Index")
forest_plot_fi_filter <- create_forest_plot(meta_results, m3sig, "fi", "the Frailty Index", T)

ggsave("forest_plot_age.png", forest_plot_age, width = 10, height = 6, units = "in")
ggsave("forest_plot_fi.png", forest_plot_fi, width = 10, height = 12, units = "in")
ggsave("forest_plot_fi_filter.png", forest_plot_fi_filter, width = 10, height = 6, units = "in")
