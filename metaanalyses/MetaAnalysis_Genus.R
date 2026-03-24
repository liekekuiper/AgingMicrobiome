rm(list = ls())
library("data.table")
library("dplyr")
library("metafor")
library("ggplot2")
library("ggrepel")
library("tidyr")
library("ComplexUpset")
library("VennDiagram")
library('stringr')
library("stringr")
library("ggpubr")


#### Data loading and preprocessing ####
# Load all datasets into a named list
files <- list(
  FHS16s = "FHS/Results_FHS_16s.csv",
  FHSshot = "FHS/Results_FHS_shot.csv",
  MrOS16s = "MrOS/Results_MrOS_16s.csv",
  MrOSshot = "MrOS/Results_MrOS_shot.csv",
  RS16s = "RS/Results_RS_16s.csv",
  DCS16s = "DCS/Results_DCS_16s.csv",
 # LLshot = "Lifelines_Frailty_results/species_qza_and_nuslurm_Script/Results_DAG3_16s/metagenomics_2025-06-07.csv",
  LLshot = "Lifelines_Frailty_results/intermediatefiles/Results_age_5_DAG3_2026-02-13.csv",
  SOL16s = "SOL/Results_SOL_16s.csv",
  SOLshot = "SOL/Results_SOL_shot.csv"
)

# Load datasets and ensure names are assigned
dataframes <- setNames(lapply(files, fread, data.table = FALSE), names(files))

names(dataframes) <- names(files)

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

# Add mean_age to each dataframe
for (name in names(dataframes)) {
  dataframes[[name]]$mean_age <- mean_cohort_age[[name]]
}


# Expand MrOS datasets for subgroup analyses
expand_MrOS <- function(df) {
  expanded <- list(
    df,
    df %>% mutate(Datasplit = "age_5"),
    df %>% mutate(Datasplit = "men")
  )
  bind_rows(expanded)
}

if ("MrOS16s" %in% names(dataframes)) {
  dataframes$MrOS16s <- expand_MrOS(dataframes$MrOS16s)
}
if ("MrOSshot" %in% names(dataframes)) {
  dataframes$MrOSshot <- expand_MrOS(dataframes$MrOSshot)
}

# Identify prefixes
prefixes <- unique(sub("(16s|shot)$", "", names(dataframes)))

# Add cohort column to each dataframe
dataframes <- mapply(function(df, name) {
  df$cohort <- sub("(16s|shot)$", "", name)
  return(df)
}, dataframes, names(dataframes), SIMPLIFY = FALSE)


# Initialize lists to hold the subsets
data_16s <- list()
data_shot <- list()
data_all <- list()

# Loop over prefixes to create `16s`, `shot`, and `all`
for (prefix in prefixes) {
  # Check if both `16s` and `shot` exist for the prefix
  df_16s <- dataframes[[paste0(prefix, "16s")]]
  df_shot <- dataframes[[paste0(prefix, "shot")]]
  
  if (!is.null(df_16s)) {
    data_16s[[paste0(prefix,"16s")]] <- df_16s %>% mutate(seq_type = "16s")
  }
  
  if (!is.null(df_shot)) {
    data_shot[[paste0(prefix,"shot")]] <- df_shot %>% mutate(seq_type = "shot")
  }
  
  # If both exist, use `shot` for `all`, otherwise use the available one
  if (!is.null(df_16s) & !is.null(df_shot)) {
    data_all[[paste0(prefix,"all")]] <- df_shot %>% mutate(seq_type = "all")
  } else if (!is.null(df_16s)) {
    data_all[[paste0(prefix,"all")]] <- df_16s %>% mutate(seq_type = "all")
  } else if (!is.null(df_shot)) {
    data_all[[paste0(prefix,"all")]] <- df_shot %>% mutate(seq_type = "all")
  }
}

# Combine results into the main dataframes list
dataframes <- c(
  data_16s = data_16s,
  data_shot = data_shot,
  data_all = data_all
)
# Remove the portion before the dot in dataframe names
names(dataframes) <- sub("^[^.]+\\.", "", names(dataframes))

# Verify the updated names
print(names(dataframes))

#### First impression harmonization ####
# Function to assign Mnumber
mutate_summ_data <- function(df) {
  df %>% mutate(
    Mnumber = case_when(
      !grepl("statin", Model) ~ "m1",
      grepl("statin", Model) & !grepl("bmi", Model) ~ "m2",
      grepl("bmi", Model) & !grepl("dietscore", Model) ~ "m3",
      grepl("dietscore", Model) ~ "m4",
      TRUE ~ NA_character_
    ) 
  ) %>% filter(., !is.na(P))
}


# Apply mutation to all relevant dataframes
dataframes <- lapply(dataframes, mutate_summ_data)

# Calculate max N for each dataframe and combine into a single dataframe
max_N_df <- bind_rows(lapply(names(dataframes), function(dataset_name) {
  df <- dataframes[[dataset_name]]  # Access the dataframe by name
  df %>%
    group_by(Mnumber, Outcome, Datasplit, seq_type) %>%
    summarise(max_N = max(N, na.rm = TRUE), .groups = "drop") %>%
    mutate(dataset = dataset_name)  # Add the dataset name for reference
}))

# Calculate the sum of max_N for each subgroup in max_N_df
sum_max_N_df <- max_N_df %>%
  group_by(Mnumber, Outcome, Datasplit, seq_type) %>%
  summarise(
    sum_max_N = sum(max_N, na.rm = TRUE), 
    nstudies = n(), 
    .groups = "drop"
  )
sum_max_N_df_all = subset(sum_max_N_df, seq_type == 'all' & Outcome %in% c('age', 'fi'))
sum_max_N_df_all_m3 = subset(sum_max_N_df_all, Mnumber == 'm3')

# Function to create a binary presence/absence dataframe for each dataset
get_risk <- function(df, risk_level) {
  df <- df %>% select(Variable) %>% distinct()  # Keep only unique Variable names
  df <- df %>%  filter(grepl("^g__", Variable))
  df[[risk_level]] <- 1  # Mark presence as 1
  return(df)
}

# Apply function to each dataset and store results
sig_risk_list <- lapply(names(dataframes), function(name) {
  df <- dataframes[[name]]
  df <- df %>% filter(Datasplit == "all")  # Apply filter before processing
  get_risk(df, name)
})


# Combine all datasets into one binary matrix
all_risk_data <- Reduce(function(x, y) full_join(x, y, by = "Variable"), sig_risk_list)

# Replace NA values with 0 (absence of variable in that dataset)
all_risk_data[is.na(all_risk_data)] <- 0

# Remove Variable column as it is not needed for the UpSet plot
rownames(all_risk_data) = all_risk_data$Variable
all_risk_data <- all_risk_data %>% select(-Variable)

# Ensure all columns are numeric (0 or 1)
all_risk_data[] <- lapply(all_risk_data, function(x) as.integer(x > 0))
all_risk_data = dplyr::select(all_risk_data, !contains("all"))
colnames(all_risk_data) = gsub('shot', ' MG', colnames(all_risk_data))
colnames(all_risk_data) = gsub('16s', ' 16S', colnames(all_risk_data))
# Define the risk level columns
risk_vars <- colnames(all_risk_data)
risk_vars_ordered <- rev(c('FHS 16S','FHS MG','MrOS 16S','MrOS MG','SOL 16S','SOL MG','DCS 16S',
                       'RS 16S','LL MG'))

# Generate the UpSet plot
upsetplotgenus = upsetplotgenus <- ComplexUpset::upset(
  all_risk_data,
  intersect = risk_vars_ordered,
  width_ratio = 0.4,
  name = "Combination",
  group_by = "degree",                      
  sort_intersections = "descending",         
  sort_intersections_by = c("degree","cardinality"),
  sort_sets = FALSE,
  matrix = (
    intersection_matrix(geom = geom_point(shape = "circle filled", size = 6.5)) +
      scale_color_manual(
        values = c(
          "FHS 16S"="#cc6677","FHS MG"="#332288","MrOS 16S"="#ddcc77","MrOS MG"="#117733",
          "RS 16S"="#88ccee","LL MG"="#882255","DCS 16S"="#44AA99","SOL 16S"="#999933","SOL MG"="#AA4499"
        ),
        guide = guide_legend(override.aes = list(shape = "circle"))
      )
  ),
  queries = list(
    upset_query(set = "FHS 16S", fill = "#cc6677"),
    upset_query(set = "FHS MG",  fill = "#332288"),
    upset_query(set = "MrOS 16S",fill = "#ddcc77"),
    upset_query(set = "MrOS MG", fill = "#117733"),
    upset_query(set = "RS 16S",  fill = "#88ccee"),
    upset_query(set = "LL MG",   fill = "#882255"),
    upset_query(set = "DCS 16S", fill = "#44AA99"),
    upset_query(set = "SOL 16S", fill = "#999933"),
    upset_query(set = "SOL MG",  fill = "#AA4499")
  ),
  themes = upset_default_themes(text = element_text(size = 45)),
  base_annotations = list("Number of core genera\noverlapping for combination" =
                            intersection_size(text = list(size = 9), counts = TRUE)),
  set_sizes = # replacing upset_set_size() to show size of sets
  upset_set_size(position = "left") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    hjust = 1.05,     # nudge labels past bar ends
    size = 9
  ) +
  ylab("Number of core\ngenera per dataset") +
  coord_cartesian(clip = "off")
)


ggsave(plot = upsetplotgenus, filename = "UpsetPlotGenus.png", width = 27, height = 15, units='in', bg = "white")
ggsave(
  plot = upsetplotgenus,
  filename = "UpsetPlotGenus.pdf",
  width = 27,
  height = 15,
  units = "in",
  bg = "white"
)

#
exclusive_DCS16s <- all_risk_data %>%
  filter(`DCS 16S` == 1 & rowSums(select(., -`DCS 16S`)) == 0)
exclusive_LLshot <- all_risk_data %>%
  filter(`LL MG` == 1 & rowSums(select(., -`LL MG`)) == 0)

shot_columns <- grep("MG$", colnames(all_risk_data), value = TRUE)
present_in_all_shot <- all_risk_data %>%
  filter(rowSums(select(., all_of(shot_columns))) == length(shot_columns))

s16_columns <- grep("16S$", colnames(all_risk_data), value = TRUE)
present_in_all_16s <- all_risk_data %>%
  filter(rowSums(select(., all_of(s16_columns))) == length(s16_columns))


#### Random effect meta analyses ####
# Initialize results list
results <- list()
mnumber_levels <- c("m1", "m2", "m3", "m4")
outcome_levels <- c("age", "fi")
population_levels <- unique(dataframes$LLshot$Datasplit)
risk_levels <- unique(unlist(lapply(dataframes, function(df) df$Variable)))
risk_levels <- risk_levels[!grepl("species|s__", risk_levels)]

# Nested loop for meta-analysis
for (outcome in outcome_levels) {
  for (population in population_levels) {
    for (mnumber in mnumber_levels) {
      for (analysis_type in c("16s", "shot", "all")) {
        for (risk in risk_levels) {
          all_studies <- lapply(dataframes, function(df) {
            if (risk %in% df$Variable) {
              df %>% filter(
                seq_type == analysis_type,
                Variable == risk,
                Mnumber == mnumber,
                Outcome == outcome,
                Datasplit == population
              ) %>% select(N, Coefficient, Std.Error, P, LL, UL)
            } else {
              NULL
            }
          }) %>% bind_rows()
          
          if (nrow(all_studies) > 1) {
            if(nrow(all_studies) > 6) {
              print("Too many studies for risk:", risk, "\noutcome:", outcome, "\nmodel:", mnumber,  "\npopulation:", population,)
              break
            }
            # Calculate ndirection
            positive_coeff <- sum(all_studies$Coefficient > 0, na.rm = TRUE)
            negative_coeff <- sum(all_studies$Coefficient < 0, na.rm = TRUE)
            ndirection <- max(positive_coeff, negative_coeff)
            
            # Try random-effects model with fallback to fixed-effects
            meta_result <- tryCatch(
              rma(yi = Coefficient, sei = Std.Error, data = all_studies, method = "REML"),
              error = function(e) {
                message(paste("Random-effects model failed for risk:", risk, "\noutcome:", outcome, "\nmodel:", mnumber,  "\npopulation:", population, '\nanalysis type:', analysis_type, ".\nSwitching to fixed-effects model."))
                rma(yi = Coefficient, sei = Std.Error, data = all_studies, method = "FE")
              }
            )            
            results <- append(results, list(
              list(
                RISK = risk,
                Mnumber = mnumber,
                OUTCOME2 = outcome,
                Population = population,
                Analysis_Type = analysis_type,
                B = meta_result$b,
                LL = meta_result$ci.lb,
                UL = meta_result$ci.ub,
                SE = meta_result$se,
                P = meta_result$pval,
                QE = meta_result$QE,
                QEp = meta_result$QEp,
                I2 = meta_result$I2,
                nstudies = nrow(all_studies),
                tau2 = meta_result$tau2,
                ndirection = ndirection,
                nparticipants = sum(all_studies$N)
              )
            ))
          }
        }
      }
    }
  }
}

# Compile results into a dataframe
meta_results_df <- do.call(rbind, results) %>%
  as.data.frame() %>%
  mutate(across(c(RISK, Mnumber, OUTCOME2, Population, Analysis_Type), as.character),
         across(c(B, LL, UL, SE, P, QE, QEp, I2, nstudies, tau2, ndirection, nparticipants), as.numeric))
#only later look up significant ones
submetacheck = subset(meta_results_df, Population != 'all' & Analysis_Type != 'all' & !OUTCOME2 %in% c('cont', "mortality"))
submetaresults = subset(meta_results_df, Population != 'all' & Analysis_Type == 'all'& !OUTCOME2 %in% c('cont', "mortality"))

meta_results_df = subset(meta_results_df, Population == 'all' & !OUTCOME2 %in% c('cont', "mortality"))
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

#m3sig16s = subset(meta_results, Mnumber == 'm3' & Analysis_Type == "16s" & Population == 'all' & pFDR < 0.05)
#m3sigshot= subset(meta_results, Mnumber == 'm3' & Analysis_Type == "shot" & Population == 'all' & pFDR < 0.05)


`%+%` <- function(a, b) paste0(a, b)

# Label cleaner functie
clean_label <- function(varname, table = FALSE) {
  if (is.na(varname)) {
    return(NA_character_)
  }
  
  if (str_starts(varname, "g__")) {
    name <- str_remove(varname, "^g__")
    genus <- case_when(
      name %in% c("CAG_217", "CAG_273") ~ "Clostridium (sp. CAG:" %+% str_extract(name, "\\d+") %+% ")",
      name == "CAG_83" ~ "Oscillibacter (sp. CAG:83)",
      name == "ER4" ~ "Oscillibacter (sp. ER4)",
      name == "SFMI01" ~ "Christensenellales (sp. SFMI01)",
      name == "UBA1417" ~ "Acutalibacter (sp. UBA1417)",
      TRUE ~ str_replace(name, "_.*$", "")
    )
    return(if (table) genus else paste0("italic('", genus, "')"))
  }
  
  get_type <- function(v) {
    case_when(
      str_detect(v, "genus") ~ "Genus",
      str_detect(v, "asv") ~ "ASV",
      str_detect(v, 'feature') ~ 'ASV',
      TRUE ~ ""
    )
  }
  
  if (str_starts(varname, "simpson_e")) {
    type <- get_type(varname)
    return(if (table) paste0("Inverse Simpson (", type, ")") else paste0("Inverse~Simpson~(", type, ")"))
  }
  
  if (str_starts(varname, "simpson")) {
    type <- get_type(varname)
    return(if (table) paste0("Simpson (", type, ")") else paste0("Simpson~(", type, ")"))
  }
  
  if (str_starts(varname, "chao1")) {
    type <- get_type(varname)
    return(if (table) paste0("Chao1 (", type, ")") else paste0("Chao1~(", type, ")"))
  }
  
  if (str_starts(varname, "shannon")) {
    type <- get_type(varname)
    return(if (table) paste0("Shannon (", type, ")") else paste0("Shannon~(", type, ")"))
  }
  
  if (str_starts(varname, "min_")) {
    beta <- case_when(
      str_detect(varname, "bray") ~ "Bray-Curtis",
      str_detect(varname, "jacc") ~ "Jaccard",
      str_detect(varname, "aitch") ~ "Aitchison",
      str_detect(varname, "uu") ~ "Unweighted UniFrac",
      str_detect(varname, "wu") ~ "Weighted UniFrac",
      TRUE ~ NA_character_
    )
    type <- get_type(varname)
    return(if (table) 
      paste0("Uniqueness based on ", beta, " (", type, ")") 
      else 
        paste0("Uniqueness~based~on~", beta, "~(", type, ")"))
  }
  
  return(varname)
}

# --------------------------------------------------
# Helper function: convert all columns to character
# and save to ../Paper/Supplements/ using fwrite
# --------------------------------------------------
write_supp_table <- function(df, filename) {
  df %>%
    mutate(across(everything(), as.character)) %>%
    data.table::fwrite(
      file = file.path("../Paper/Supplements", filename),
      na = "",
      quote = TRUE
    )
}

supptabl2 = m3_meta_results_df %>% subset(., OUTCOME2 == 'age') %>%
  arrange(pFDR) %>% 
  mutate(
    Variable = sapply(RISK, clean_label, table = TRUE),
    # Confidence intervals for original values
    B_CI = case_when(OUTCOME2 == 'fi' ~ paste0(sprintf(B, fmt = "%.4f"), " (", sprintf(LL, fmt = "%.4f"), "; ", sprintf(UL, fmt = "%.4f"), ")"),
                              TRUE ~ paste0(sprintf(B, fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), "; ", sprintf(UL, fmt = "%.2f"), ")")),
    # Format p-values
    pFDR_handig = case_when(pFDR < 0.01 ~ sub("e","x10^",sprintf(pFDR, fmt="%.2e")),
                            pFDR >= 0.01 ~ sprintf(pFDR, fmt = "%.2f")),
    QE = sprintf(QE, fmt = "%.2f"),
    QEp_handig = case_when(QEp < 0.01 ~ sub("e","x10^",sprintf(QEp, fmt="%.2e")),
                           QEp >= 0.01 ~ sprintf(QEp, fmt = "%.2f")),
    I2 = sprintf(I2, fmt= "%.2f"),
    tau2 = sprintf(tau2, fmt= "%.2f"),
    Outcome = case_when(OUTCOME2 == 'age' ~ 'Chronological age',
                        OUTCOME2 == 'fi' ~ 'Frailty index'),
    Studiesdir = paste0(nstudies, " (", ndirection, ")")
  )  %>%
  select(Outcome, Variable, RISK, Studiesdir, nparticipants, B_CI, pFDR_handig, QE, QEp_handig, I2, tau2) %>%
  rename(., 
         `Microbial Feature` = Variable, 
         `GreenGenes2 name` = RISK,
         `Nstudies (Ndirection)`= Studiesdir,
         Nparticipants = nparticipants,
         `Beta (CI)` = B_CI,
         pFDR = pFDR_handig,
         `QE p-value` =QEp_handig)
write_supp_table(supptabl2, "SuppTable2.csv")

supptabl8 = m3_meta_results_df %>% subset(., OUTCOME2 == 'fi') %>%
  arrange(pFDR) %>% 
  mutate(
    Variable = sapply(RISK, clean_label, table = TRUE),
    # Confidence intervals for original values
    B_CI = case_when(OUTCOME2 == 'fi' ~ paste0(sprintf(B, fmt = "%.4f"), " (", sprintf(LL, fmt = "%.4f"), "; ", sprintf(UL, fmt = "%.4f"), ")"),
                     TRUE ~ paste0(sprintf(B, fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), "; ", sprintf(UL, fmt = "%.2f"), ")")),
    # Format p-values
    pFDR_handig = case_when(pFDR < 0.01 ~ sub("e","x10^",sprintf(pFDR, fmt="%.2e")),
                            pFDR >= 0.01 ~ sprintf(pFDR, fmt = "%.2f")),
    QE = sprintf(QE, fmt = "%.2f"),
    QEp_handig = case_when(QEp < 0.01 ~ sub("e","x10^",sprintf(QEp, fmt="%.2e")),
                           QEp >= 0.01 ~ sprintf(QEp, fmt = "%.2f")),
    I2 = sprintf(I2, fmt= "%.2f"),
    tau2 = sprintf(tau2, fmt= "%.2f"),
    Outcome = case_when(OUTCOME2 == 'age' ~ 'Chronological age',
                        OUTCOME2 == 'fi' ~ 'Frailty index'),
    Studiesdir = paste0(nstudies, " (", ndirection, ")")
  )  %>%
  select(Outcome, Variable, RISK, Studiesdir, nparticipants, B_CI, pFDR_handig, QE, QEp_handig, I2, tau2) %>%
  rename(.,
         `Microbial Feature` = Variable, 
         `GreenGenes2 name` = RISK,
         `Nstudies (Ndirection)`= Studiesdir,
         Nparticipants = nparticipants,
         `Beta (CI)` = B_CI,
         pFDR = pFDR_handig,
         `QE p-value` =QEp_handig)
write_supp_table(supptabl8, "SuppTable8.csv")

dietsens = meta_results %>%
  subset(Mnumber == 'm4' & Analysis_Type == 'all' & paste(RISK, OUTCOME2) %in% paste(m3sig$RISK, m3sig$OUTCOME2)) %>%
  mutate(
    Variable = sapply(RISK, clean_label, table = TRUE),
    # Confidence intervals for original values
    B_CI = case_when(OUTCOME2 == 'fi' ~ paste0(sprintf(B, fmt = "%.4f"), " (", sprintf(LL, fmt = "%.4f"), "; ", sprintf(UL, fmt = "%.4f"), ")"),
                     TRUE ~ paste0(sprintf(B, fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), "; ", sprintf(UL, fmt = "%.2f"), ")")),
    # Format p-values
    pFDR_handig = case_when(pFDR < 0.01 ~ sub("e","x10^",sprintf(pFDR, fmt="%.2e")),
                            pFDR >= 0.01 ~ sprintf(pFDR, fmt = "%.2f")),
    QE = sprintf(QE, fmt = "%.2f"),
    QEp_handig = case_when(QEp < 0.01 ~ sub("e","x10^",sprintf(QEp, fmt="%.2e")),
                           QEp >= 0.01 ~ sprintf(QEp, fmt = "%.2f")),
    I2 = sprintf(I2, fmt= "%.2f"),
    tau2 = sprintf(tau2, fmt= "%.2f"),
    
    Studiesdir = paste0(nstudies, " (", ndirection, ")"),
    Outcome = case_when(OUTCOME2 == 'age' ~ 'Chronological age',
                        OUTCOME2 == 'fi' ~ 'Frailty index')
  )  %>%
  select(Outcome, Variable, RISK, Studiesdir, nparticipants, B_CI, pFDR_handig, QE, QEp_handig, I2, tau2) %>%
  rename(., 
         `Microbial Feature` = Variable, 
         `GreenGenes2 name` = RISK,
         `Nstudies (Ndirection)`= Studiesdir,
         Nparticipants = nparticipants,
         `Beta (CI)` = B_CI,
         pFDR = pFDR_handig,
         `QE p-value` =QEp_handig)
write_supp_table(dietsens, "SuppTable3.csv")

m3display = m3sig %>%
  mutate(
    # Store original values before transformation
    B_original = B,
    LL_original = LL,
    UL_original = UL,
    
    # Apply scaling for B, LL, and UL
    across(c(B, LL, UL), ~ case_when(
      grepl("min_", RISK) & OUTCOME2 == "age" ~ .x * 0.1,
      grepl("min_", RISK) & OUTCOME2 == "fi" ~ .x * 10,
      OUTCOME2 == "fi" ~ .x * 100,
      TRUE ~ .x
    )),
    Variable = sapply(RISK, clean_label, table = TRUE),
    # Confidence intervals for scaled values
    B_CI = paste0(sprintf(B, fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), "; ", sprintf(UL, fmt = "%.2f"), ")"),
    
    # Confidence intervals for original values
    B_CI_original = case_when(OUTCOME2 == 'fi' ~ paste0(sprintf(B_original, fmt = "%.4f"), " (", sprintf(LL_original, fmt = "%.4f"), "; ", sprintf(UL_original, fmt = "%.4f"), ")"),
                              TRUE ~ paste0(sprintf(B_original, fmt = "%.2f"), " (", sprintf(LL_original, fmt = "%.2f"), "; ", sprintf(UL_original, fmt = "%.2f"), ")")),
    
    # Format p-values
    pFDR_handig = case_when(pFDR < 0.01 ~ sub("e","x10^",sprintf(pFDR, fmt="%.2e")),
                            pFDR >= 0.01 ~ sprintf(pFDR, fmt = "%.2f")),
    
    QE = sprintf(QE, fmt = "%.2f"),
    
    QEp_handig = case_when(QEp < 0.01 ~ sub("e","x10^",sprintf(QEp, fmt="%.2e")),
                           QEp >= 0.01 ~ sprintf(QEp, fmt = "%.2f")),
    
    I2 = sprintf(I2, fmt= "%.2f"),
    tau2 = sprintf(tau2, fmt= "%.2f"),
    
    Studiesdir = paste0(nstudies, " (", ndirection, ")")
  ) %>%
  select(OUTCOME2, Variable, RISK, Studiesdir, nparticipants,B_original, B_CI, B_CI_original, pFDR_handig, QE, QEp_handig, I2, tau2, pFDR)



# Save feature-genus mapping
val_genus = select(m3display, c(RISK, OUTCOME2)) %>%
  rename(., 'Feature_Genus'='RISK', 'Outcome'='OUTCOME2')
data.table::fwrite(val_genus, "Updated_CHARGE_genus.csv")

# Create separate displays for age and fi outcomes with both original and scaled CIs
m3agedisplay = m3display %>% 
  filter(OUTCOME2 == 'age') %>%
  select(Variable,RISK, Studiesdir, nparticipants, B_CI, B_CI_original, pFDR_handig, QE, QEp_handig, I2, tau2, pFDR)

m3betafidisplay = m3display %>% 
  filter(OUTCOME2 == 'fi') %>%
  select(Variable,RISK, Studiesdir, nparticipants, B_original, B_CI, B_CI_original, pFDR_handig, QE, QEp_handig, I2, tau2, pFDR) %>%
  rename(., B = B_original)
m3fidisplay = select(m3betafidisplay, -B)

# Save updated files
write.csv(m3agedisplay, "M3AgeDisplay.csv", row.names = FALSE)
write.csv(m3fidisplay, "M3FIDisplay.csv", row.names = FALSE)
write.csv(m3betafidisplay, "M3_B_FIDisplay.csv", row.names = FALSE)
#### Subgroup analyses ####
#### Heterogeneity
library(patchwork)

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

# Subset heterogeneous associations
heterogeneous <- subset(m3sig, QEp < 0.05 | I2 >= 40)

# Unique combinations of RISK and OUTCOME2
heterogeneous_combinations <- unique(heterogeneous[, c("RISK", "OUTCOME2")])

# Define cohorts used in meta-analyses
dataframes_all <- list(FHS = dataframes$FHSall, 
                       MrOS = dataframes$MrOSall, 
                       RS = dataframes$RSall, 
                       DCS = dataframes$DCSall, 
                       LL = dataframes$LLall, 
                       SOL = dataframes$SOLall)

# Function to extract relevant data for each cohort
extract_cohort_data <- function(df, cohort_name, risk, outcome) {
  df_sub <- df %>%
    filter(Mnumber == "m3", Variable == risk, Outcome == outcome, Datasplit == 'all')
  
  if (nrow(df_sub) == 0) return(NULL)
  
  df_sub %>%
    mutate(Cohort = paste0(cohort_name, " (N=", N, ")"),
           B = Coefficient) %>%
    select(Cohort, B, LL, UL)
}

# Generate individual forest plots
forest_plot_list <- list()


for (i in 1:nrow(heterogeneous_combinations)) {
  
  risk <- heterogeneous_combinations$RISK[i]
  outcome <- heterogeneous_combinations$OUTCOME2[i]
  outcome2 <- ifelse(outcome == 'age', 'Chronological Age',
                     ifelse(outcome == 'fi', 'the Frailty Index', outcome))
  
  # Meta-analysis line
  meta_effect <- heterogeneous %>%
    filter(RISK == risk, OUTCOME2 == outcome) %>%
    mutate(Cohort = paste0('Meta-analysis', " (N=", nparticipants, ")")) %>%
    select(Cohort, B, LL, UL)
  
  # Extract per-cohort results
  cohort_results <- do.call(
    rbind,
    Filter(Negate(is.null),
           lapply(names(dataframes_all), function(cohort) {
             extract_cohort_data(dataframes_all[[cohort]], cohort, risk, outcome)
           })
    )
  )
  
  if (nrow(cohort_results) == 0) next
  
  forest_data <- rbind(cohort_results, meta_effect)
  meta_name <- unique(meta_effect$Cohort)
  
  # Order full labels for display
  forest_data$Cohort <- factor(
    forest_data$Cohort,
    levels = c(setdiff(unique(forest_data$Cohort), meta_name), meta_name)
  )

  # Add short names for consitent colors
  forest_data <- forest_data %>%
    mutate(
      CohortShort = ifelse(
        grepl("^Meta-analysis", Cohort),
        "Meta-analysis",
        sub(" \\(N=.*\\)$", "", Cohort)
      )
    )
  
  present_levels <- intersect(cohort_names, unique(forest_data$CohortShort))
  present_levels <- c(
    setdiff(present_levels, "Meta-analysis"),
    intersect("Meta-analysis", present_levels)
  )
  
  forest_data$CohortShort <- factor(forest_data$CohortShort, levels = present_levels)
  
  # Clean label
  label_expr <- tryCatch({
    parse(text = clean_label(risk, table = FALSE))[[1]]
  }, error = function(e) { risk })
  
  # Build forest plot

  plot <- ggplot(forest_data, aes(y = Cohort, x = B, xmin = LL, xmax = UL)) +
    geom_pointrange(aes(color = CohortShort), size = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    scale_color_manual(values = cohort_colours) +
    theme_minimal(base_size = 15) +
    labs(
      title = label_expr,
      x = "Effect Size (B)",
      y = "Cohort"
    ) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(
        face = ifelse(levels(forest_data$Cohort) == meta_name, "bold", "plain")
      )
    )
  
  forest_plot_list[[paste(risk, outcome, sep = "_")]] <- plot
}

#Save individual plots
for (name in names(forest_plot_list)) {
  ggsave(
    paste0(name, "_forestplot.png"),
    forest_plot_list[[name]],
    width = 8, height = 6, units = "in", limitsize = FALSE
  )
}

#Group plots by outcome
plots_by_outcome <- list()

for (name in names(forest_plot_list)) {
  suffix <- sub(".*_", "", name)
  
  if (!suffix %in% names(plots_by_outcome)) {
    plots_by_outcome[[suffix]] <- list()
  }
  
  if (inherits(forest_plot_list[[name]], "gg")) {
    plots_by_outcome[[suffix]] <- append(
      plots_by_outcome[[suffix]],
      list(forest_plot_list[[name]])
    )
  }
}


# Patchwork combinations
for (outcome in names(plots_by_outcome)) {
  
  if (length(plots_by_outcome[[outcome]]) == 0) next
  
  combined_plot <- wrap_plots(
    plots_by_outcome[[outcome]],
    ncol = 2,
    tag_level = 'new'
  ) +
    plot_annotation(tag_levels = 'a')
  
  ggsave(
    paste0("patchwork_", outcome, ".png"),
    combined_plot,
    width = 15,
    height = 5 * ceiling(length(plots_by_outcome[[outcome]]) / 2),
    units = "in"
  )
  
  ggsave(
    paste0("patchwork_", outcome, ".pdf"),
    combined_plot,
    width = 15,
    height = 5 * ceiling(length(plots_by_outcome[[outcome]]) / 2),
    units = "in"
  )
}


# Unique combinations of RISK and OUTCOME2
risk_outcome_combinations <- unique(m3sig[, c("RISK", "OUTCOME2")])


subgroupresults <- list()

# Meta-regression loop for each RISK-OUTCOME2 combination
for (i in 1:nrow(risk_outcome_combinations)) {
  risk <- risk_outcome_combinations$RISK[i]
  outcome <- risk_outcome_combinations$OUTCOME2[i]
  
  for (mnumber in "m3") {
    
    # Grouping by Datasplit categories: "age_" and "men/women"
    for (population_group in list(c("age_3", "age_2", "age_1", "age_4", "age_5"), c("women", "men"))) {
      population_filter <- unlist(population_group)
      
      # Filter relevant studies based on the current group
      all_studies <- lapply(dataframes_all, function(df) {
        df_filtered <- df %>%
          filter(Mnumber == mnumber,
                 Outcome == outcome,
                 Datasplit %in% population_filter)
        
        if (risk %in% df_filtered$Variable) {
          df_filtered <- df_filtered %>% filter(Variable == risk)
          return(df_filtered %>% select(N, Coefficient, Std.Error, P, LL, UL, Datasplit, cohort))
        } else {
          return(NULL)
        }
      }) %>% bind_rows()
      
      # Run meta-regression if data is available
      if (nrow(all_studies) > 1) {
        
        # Convert Datasplit to a factor and set reference level
        all_studies <- all_studies %>%
          mutate(Datasplit = factor(Datasplit, levels = if ("age_3" %in% population_filter) {
            c("age_3", setdiff(population_filter, "age_3"))
          } else {
            c("women", setdiff(population_filter, "women"))
          }))
        ref_level <- levels(all_studies$Datasplit)[1]
        
        # Perform meta-regression with Datasplit as moderator
        meta_result <- tryCatch(
          rma(yi = Coefficient, sei = Std.Error, mods = ~ Datasplit, data = all_studies, method = "REML"),
          error = function(e) {
            message(paste("Random-effects model failed for risk:", risk, ". Switching to fixed-effects model."))
            rma(yi = Coefficient, sei = Std.Error, mods = ~ Datasplit, data = all_studies, method = "FE")
          }
        )
        print(risk_outcome_combinations[i,])
        print(meta_result)
      
        # Extract coefficients for each subgroup
        estimates <- as.data.frame(meta_result$b)
        colnames(estimates) <- c("B")
        estimates$SE <- meta_result$se
        estimates$P <- meta_result$pval
        estimates$LL <- meta_result$ci.lb
        estimates$UL <- meta_result$ci.ub
        estimates$Datasplit <- rownames(estimates)
        
        # Store each subgroup separately in a new row
        for (j in 1:nrow(estimates)) {
          subgroup_name <- gsub("Datasplit", "", estimates$Datasplit[j])
          if (subgroup_name == "intrcpt") {
            subgroup_name <- ref_level
            # If it's the intercept, include all unique cohorts for nstudies
            nstudies_subgroup <- length(unique(all_studies$cohort))
            nparticipants_subgroup <- sum(all_studies$N, na.rm = TRUE)
            ndirection_subgroup = '-'
          } else {
            # Otherwise, calculate per Datasplit subgroup
            nstudies_subgroup <- nrow(all_studies %>% filter(Datasplit == subgroup_name))
            nparticipants_subgroup <- all_studies %>%
              filter(Datasplit == subgroup_name) %>%
              summarise(nparticipants = sum(N, na.rm = TRUE)) %>%
              pull(nparticipants)
            subgroup_data <- all_studies %>% filter(Datasplit == subgroup_name)
            positive_coeff <- sum(subgroup_data$Coefficient > 0, na.rm = TRUE)
            negative_coeff <- sum(subgroup_data$Coefficient < 0, na.rm = TRUE)
            ndirection_subgroup <- max(positive_coeff, negative_coeff)
          }
          
          subgroupresults <- append(subgroupresults, list(
            list(
              RISK = risk,
              OUTCOME2 = outcome,
              Mnumber = mnumber,
              Population_Group = subgroup_name,
              B = estimates$B[j],
              SE = estimates$SE[j],
              P = estimates$P[j],
              LL = estimates$LL[j],
              UL = estimates$UL[j],
              QE = meta_result$QE,
              QEp = meta_result$QEp,
              I2 = meta_result$I2,
              nstudies = nstudies_subgroup,
              tau2 = meta_result$tau2,
              nparticipants = nparticipants_subgroup,
              QM = meta_result$QM,
              QMp = meta_result$QMp,
              ndirection = ndirection_subgroup,
              R2 = meta_result$R2
            )
          ))
        }
      }
    }
  }
}

# Convert list to dataframe
meta_regression_df <- do.call(rbind, subgroupresults) %>%
  as.data.frame() %>%
  mutate(
    P = as.numeric(P),
    QMp = as.numeric(QMp)
  )

#Every test is now in here several times, group them so we only FDR correct for number of test we have run
unique_tests <- meta_regression_df %>%
  select(RISK, OUTCOME2, Mnumber, I2, tau2, QM, QE, QMp, QEp) %>%
  distinct() %>%
  mutate(
    QMpFDR = p.adjust(QMp, method = "fdr"),
    QEpFDR = p.adjust(QEp, method = "fdr")
  )


meta_fdr <- left_join(meta_regression_df, unique_tests_fdr) %>%
  mutate(., pFDR = p.adjust(P, method = 'fdr'))

meta_fdr_clean <- meta_fdr %>%
  group_by(RISK, OUTCOME2, Mnumber, I2, tau2) %>%
  mutate(
    # Identify a single row that will display QM/QMpFDR/QEp/QEpFDR
    show_tests = row_number() == 2,
    
    QM_display    = ifelse(show_tests, QM, NA),
    QMpFDR_display = ifelse(show_tests, QMpFDR, NA),
    QE_display     = ifelse(show_tests, QE, NA),
    QEpFDR_display = ifelse(show_tests, QEpFDR, NA),
    R2_display = case_when(QMpFDR < 0.05 & show_tests ~ sprintf("%.2f", as.numeric(R2)),
                           QMpFDR >= 0.05 & show_tests ~ "-",
                           !show_tests ~ NA_character_)
    
  ) %>%
  ungroup()

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
    pFDR_handig = case_when(pFDR < 0.01 ~ sub("e","x10^",sprintf(pFDR, fmt="%.2e")),
                            pFDR >= 0.01 ~ sprintf(pFDR, fmt = "%.2f")),
    QM_display = case_when(is.na(QM_display) ~ NA_character_,
                           TRUE ~ sprintf(QM_display, fmt = "%.2f")),
    QE_display = case_when(is.na(QE_display) ~ NA_character_, 
                           TRUE ~ sprintf(QE_display, fmt = "%.2f")),
    QEp_handig = case_when(QEpFDR_display < 0.01 ~ sub("e","x10^",sprintf(QEpFDR_display, fmt="%.2e")),
                           QEpFDR_display >= 0.01 ~ sprintf(QEpFDR_display, fmt = "%.2f")),
    QMpFDR_handig = case_when(QMpFDR_display < 0.01 ~ sub("e","x10^",sprintf(QMpFDR_display, fmt="%.2e")),
                              QMpFDR_display >= 0.01 ~ sprintf(QMpFDR_display, fmt = "%.2f")),
    I2 = sprintf(I2, fmt= "%.2f"),
    tau2 = sprintf(tau2, fmt= "%.2f"),
    Outcome = case_when(OUTCOME2 == 'age' ~ 'Chronological age',
                        OUTCOME2 == 'fi' ~ 'Frailty index'),
    Studiesdir = paste0(nstudies, " (", ndirection, ")")
  )  %>%
  select(Outcome, Comparison, Variable, RISK, Studiesdir, nparticipants, B_CI, pFDR_handig, QM_display, QMpFDR_handig, R2_display, QE_display, QEp_handig, I2, tau2) %>%
  rename(., 
         `Microbial Feature` = Variable,
         `GreenGenes2 name` = RISK,
         `Nstudies (Ndirection)`= Studiesdir,
         Nparticipants = nparticipants,
         `Beta (CI)` = B_CI,
         pFDR = pFDR_handig,
         QM = QM_display,
        `QM pFDR` = QMpFDR_handig,
        R2 = R2_display,
        QE = QE_display,
         `QE p-value` =QEp_handig) 
write_supp_table(supptable4, "SuppTable4.csv")

#### Figures ####
# Preprocess
process_data <- function(data, outcome, model) {
  data %>%
    mutate(
      diffexpressed = case_when(
        B > 0 & pFDR < 0.05 & ndirection > 1 ~ "UP",
        B < 0 & pFDR < 0.05 & ndirection > 1 ~ "DOWN",
        TRUE ~ "NO"
      ),
      delabel = ifelse(diffexpressed != "NO", RISK, NA)
    ) %>%
    filter(OUTCOME2 == outcome, Mnumber == model)
}

# Subsets
meta_results_all <- meta_results %>% filter(Analysis_Type == 'all')
metaage <- process_data(meta_results_all, "age", "m3")
metafi  <- process_data(meta_results_all, "fi",  "m3")
metaagem1 <- process_data(meta_results_all, "age", "m1") %>% filter(Population == "all")
metafim1  <- process_data(meta_results_all, "fi",  "m1") %>% filter(Population == "all")
meta_results_age_m3 <- filter(metaage, Population == "all")
meta_results_fi_m3 <- filter(metafi, Population == "all")

validated_genus <- fread("Validationgenus.csv") %>%
  rename(val_pFDR = pFDR,
         val_abprev = abprev)
validated_genus_fi = subset(validated_genus, Outcome == 'mortality')

# load and summarise validation table safely
validated_genus_raw <- fread("Validationgenus.csv")

validation_summary_age <- subset(validated_genus_raw, Outcome == 'age') %>%
  # make sure pFDR is numeric and abprev is character
  mutate(pFDR = as.numeric(pFDR), abprev = as.character(abprev)) %>%
  group_by(Variable) %>%
  summarize(
    val_n_valid = n(),                              # how many validation rows (sets)
    val_n_sig = sum(pFDR < 0.05, na.rm = TRUE),     # how many times pFDR < 0.05
    val_any_sig = val_n_sig > 0,                    # any significant
    val_any_lenient_sig = any(pFDR < 0.05 & abprev == "lenient", na.rm = TRUE),
    .groups = "drop"
  )

validation_summary_fi <- subset(validated_genus_raw, Outcome == 'mortality') %>%
  # make sure pFDR is numeric and abprev is character
  mutate(pFDR = as.numeric(pFDR), abprev = as.character(abprev)) %>%
  group_by(Variable) %>%
  summarize(
    val_n_valid = n(),                              # how many validation rows (sets)
    val_n_sig = sum(pFDR < 0.05, na.rm = TRUE),     # how many times pFDR < 0.05
    val_any_sig = val_n_sig > 0,                    # any significant
    val_any_lenient_sig = any(pFDR < 0.05 & abprev == "lenient", na.rm = TRUE),
    .groups = "drop"
  )

create_volcano_plot <- function(data, title, validated_genus_df) {
  # join validation summary by RISK -> Variable
  data <- data %>%
    left_join(validated_genus_df, by = c("RISK" = "Variable")) %>%
    mutate(
      # base label from your clean_label() -> plotmath-ready string
      label_base = ifelse(!is.na(RISK), sapply(RISK, clean_label), NA_character_),
      
      # Construct delabel using the new summarized rules and priorities:
      # 1) any lenient significant -> dagger (†)
      # 2) appears multiple times and all occurrences significant -> *
      # 3) any significant (but not 1 or 2) -> dubble dagger (‡)
      delabel = case_when(
        diffexpressed == "NO" ~ NA_character_,
        
        !is.na(val_any_lenient_sig) & val_any_lenient_sig & !is.na(label_base) ~
          paste0(label_base, "^bold('‡')"),
        
        !is.na(val_n_valid) & val_n_valid >= 2 & !is.na(val_n_sig) & val_n_sig == val_n_valid & !is.na(label_base) ~
          paste0(label_base, "^bold('*')"),
        
        !is.na(val_n_sig) & val_n_sig == 1 & !is.na(label_base) ~
          paste0(label_base, "^bold('†')"),
        
        TRUE ~ NA_character_
      )
      
    )
  
  # X-axis limit determination (unchanged)
  filtered_data <- data %>% filter(diffexpressed != "NO")
  if (nrow(filtered_data) > 0) {
    extreme_B <- filtered_data$B[which.max(abs(filtered_data$B))]
    max_abs_x <- abs(extreme_B) + 0.1
  } else {
    max_abs_x <- max(abs(data$B), na.rm = TRUE)
  }
  
  # Plot (parse = TRUE so the plotmath strings are interpreted)
  ggplot(data, aes(x = B, y = -log10(P), col = diffexpressed)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_text_repel(
      data = subset(data, diffexpressed != "NO" & !is.na(delabel)),
      aes(label = delabel),
      parse = TRUE,
      size = 5,
      box.padding = 0.6,
      point.padding = 0.6,
      max.overlaps = 200
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "black", "UP" = "red")) +
    labs(title = title, x = "Beta", y = "-Log10(P-value)") +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    coord_cartesian(xlim = c(-max_abs_x, max_abs_x))
}



# Create plots
plots <- list(
  Age = create_volcano_plot(meta_results_age_m3, "Age", validation_summary_age),
  FI = create_volcano_plot(meta_results_fi_m3, "Frailty Index", validation_summary_fi),
  Agem1 = create_volcano_plot(metaagem1, "Age M1", validation_summary_age),
  FIm1 = create_volcano_plot(metafim1, "Frailty Index M1", validation_summary_fi)
)

# Save plots
lapply(names(plots), function(name) {
  ggsave(paste0(name, "plot.png"), plots[[name]], width = 9, height = 5.5, units = "in")
})

# Combine plot
plotpaper <- ggarrange(plots[['Age']], plots[['FI']], nrow = 2, labels = c('a', 'b'), heights = c(1,1.6))
ggsave('CombinedPlot.png', plotpaper, width = 14, height = 12, units = "in", bg = "white")

### Subgroup
# Filter relevant studies based on the current group
plot_microbiome_effect <- function(variable, outcome, datasplit_groups, filename, line = TRUE) {
  # Combine all cohorts into a single dataframe
  df_results <- bind_rows(dataframes_all) %>%
    filter(Mnumber == 'm3', Outcome == outcome, Datasplit %in% datasplit_groups)
  
  # Filter only for the selected variable
  if (variable %in% df_results$Variable) {
    df_results <- df_results %>% filter(Variable == variable)
  } else {
    return(NULL)
  }
  
  # Assign numeric values to age categories
  df_results <- df_results %>%
    mutate(
      age_use = case_when(
        !grepl("age_", Datasplit) ~ mean_age,  # Use mean age if not in an age_ category
        Datasplit == "age_1" ~ 30,
        Datasplit == "age_2" ~ 45,
        Datasplit == "age_3" ~ 55,
        Datasplit == "age_4" ~ 65,
        Datasplit == "age_5" ~ 75
      )
    )
  
  # Clean label for title as expression
  title_label <- parse(text = clean_label(variable, table = FALSE))
  
  # Base plot with points colored by cohort
  baseplot <- ggplot(df_results, aes(x = age_use, y = Coefficient, size = N, color = cohort)) +
    geom_point(alpha = 0.7) +
    scale_x_continuous(breaks = c(30, 45, 55, 65, 75), labels = c("<40", "40–50", "50–60", "60–70", "70<")) +
    geom_hline(yintercept = 0, color = "grey60") + 
    xlab("Age category") + 
    ylab(paste(str_to_title(outcome), "coefficient")) +
    labs(
      colour = "Cohort",
      size = "N participants"
    ) +
    ggtitle(label = title_label)+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(size = guide_legend(override.aes = list(alpha = 1))) +
    NULL
  
  # Add line if requested
  if (line) {
    if (grepl("age_", datasplit_groups[1])) {
      plot <- baseplot +
        stat_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE, size = 1, color = "black")
    } else {
      plot <- baseplot +
        stat_smooth(aes(group = Datasplit, color = Datasplit),
                    method = "lm", formula = y ~ x, se = TRUE, size = 1)
    }
  } else {
    plot <- baseplot
  }
  
  # Save the plot
  ggsave(filename, plot, dpi = 1200, width = 7, height = 5, units = "in", limitsize = FALSE)
  
  return(list(plot = plot, results = df_results))
}




faecalibacterium_age_age = plot_microbiome_effect(
  variable = "g__Faecalibacterium",
  outcome = "age",
  datasplit_groups = c("age_3", "age_2", "age_1", "age_4", "age_5"),
  filename = "Faecalibacterium_agegroups.png"
)


## Make figure sequencing type

# Function to create a forest plot for a given outcome (e.g., "age" or "fi")
create_forest_plot <- function(outcome, display_name, filter=F) {
  # Filter meta_results for the specified outcome and Population == "all"
  filtered_meta <- meta_results %>%
    filter(OUTCOME2 == outcome, Population == "all", Mnumber == 'm3') 
  
  # Merge with m3sig to keep only significant risks
  filtered_meta <- filtered_meta %>%
    semi_join(m3sig[, c("RISK", "OUTCOME2")], by = c("RISK", "OUTCOME2")) %>%
    mutate(B = ifelse(grepl("min_|shannon|simpson|choa1", RISK), B * 0.1, B),
           LL = ifelse(grepl("min_|shannon|simpson|choa1", RISK), LL * 0.1, LL),
           UL = ifelse(grepl("min_|shannon|simpson|choa1", RISK), UL * 0.1, UL),
           Analysis_Type = factor(Analysis_Type, levels = c("shot", "16s", "all")))
  if(filter){
    het = heterogeneous %>%
      filter(OUTCOME2 == outcome)
    filtered_meta <- filtered_meta %>%
      filter(RISK %in% het$RISK) 
  }
  
  
  # Order by RISK first, then Analysis_Type
  filtered_meta <- filtered_meta %>%
    arrange(RISK, Analysis_Type) %>%
    mutate(y_position = row_number())  # Unique row number for y-axis
  
  # Label the row belonging to 'all'
  filtered_meta <- filtered_meta %>%
    arrange(RISK, Analysis_Type) %>%
    mutate(y_position = row_number()) %>%
    mutate(y_label = ifelse(
      Analysis_Type == "all",
      sapply(RISK, clean_label, table = FALSE),
      NA_character_
    ))
  
  
  
  # Create forest plot
  forest_plot <- ggplot(filtered_meta, aes(x = B, y = factor(y_position), xmin = LL, xmax = UL, color = Analysis_Type)) +
    geom_point(size = 3) +
    geom_errorbarh(height = 0.2) +
    labs(x = "Effect Estimate", y = NULL, title = paste("Forest Plot for", display_name)) +
    geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', alpha = 0.5) +
    theme_minimal() +
    scale_y_discrete(labels = function(x) parse(text = filtered_meta$y_label[match(x, filtered_meta$y_position)])) +
    scale_color_manual(
      name = "Included datasets",  # Rename legend title
      values = c("shot" = "#648fff", "16s" = "#ffb000", "all" = "#dc267f"),  # Customize colors
      labels = c("shot" = "Metagenomics", "16s" = "16s", "all" = "All (including metagenomics\n if available)"),
      breaks = rev(levels(filtered_meta$Analysis_Type))  # Maintain order in the legend as in the plot
    ) +
    theme(axis.text.y = element_text(size = 12),  
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          text = element_text(size = 16)) 
  
  
  return(forest_plot)
}

# Generate and save plots for 'age' and 'fi'
forest_plot_age <- create_forest_plot("age", "Chronological age")
forest_plot_fi <- create_forest_plot("fi", "the Frailty Index")
forest_plot_fi_filter <- create_forest_plot("fi", "the Frailty Index", T)

ggsave("forest_plot_age.png", forest_plot_age, width = 10, height = 6, units = "in")
ggsave("forest_plot_fi.png", forest_plot_fi, width = 10, height = 12, units = "in")
ggsave("forest_plot_fi_filter.png", forest_plot_fi_filter, width = 10, height = 6, units = "in")

