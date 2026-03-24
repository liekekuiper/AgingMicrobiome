rm(list = ls())
library("data.table")
library("dplyr")
library("metafor")
library("ggplot2")
library("ggrepel")

# Load all datasets into a named list
files <- list(
  FHSshot = "FHS/Results_species_FHS_shot.csv",
  MrOSshot = "MrOS/Results_species_MrOS_shot.csv",
  LLshot = "Lifelines_Frailty_results/Results_species_LL.csv",
  SOLshot = "SOL/Results_species_SOL_shot.csv"
)

# Load datasets
dataframes <- setNames(lapply(files, fread, data.table = FALSE), names(files))

# Duplicated rows in LL remove
#dataframes$LLshot <- dataframes$LLshot[!duplicated(dataframes$LLshot), ]

# Expand MrOS datasets for subgroup analyses
expand_MrOS <- function(df) {
  expanded <- list(
    df,
    df %>% mutate(Datasplit = "age_5"),
    df %>% mutate(Datasplit = "men")
  )
  bind_rows(expanded)
}
if ("MrOSshot" %in% names(dataframes)) {
  dataframes$MrOSshot <- expand_MrOS(dataframes$MrOSshot)
}

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
  ) %>% filter(!is.na(P))
}
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

# Meta-analysis setup
results <- list()
mnumber_levels <- c("m1", "m2", "m3", "m4")
outcome_levels <- c("age", "fi")
population_levels <- "all"

# Nested loop for meta-analysis
for (outcome in outcome_levels) {
  risk_levels <- if (outcome == "age") risk_levels_age else risk_levels_fi
  
  for (population in population_levels) {
    for (mnumber in mnumber_levels) {
      for (risk in risk_levels) {
        all_studies <- lapply(dataframes, function(df) {
          if (risk %in% df$Variable) {
            df %>% filter(
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
          if(nrow(all_studies) > 4){
            print(all_studies)
          }
          positive_coeff <- sum(all_studies$Coefficient > 0, na.rm = TRUE)
          negative_coeff <- sum(all_studies$Coefficient < 0, na.rm = TRUE)
          ndirection <- max(positive_coeff, negative_coeff)
          
          meta_result <- tryCatch(
            rma(yi = Coefficient, sei = Std.Error, data = all_studies, method = "REML"),
            error = function(e) {
              message(paste("Random-effects model failed for risk:", risk, ". Switching to fixed-effects model."))
              rma(yi = Coefficient, sei = Std.Error, data = all_studies, method = "FE")
            }
          )
          
          results <- append(results, list(
            list(
              RISK = risk,
              Mnumber = mnumber,
              OUTCOME2 = outcome,
              Population = population,
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

# Compile results into a dataframe and adjust p-values
meta_results_df <- do.call(rbind, results) %>%
  as.data.frame() %>%
  mutate(across(c(RISK, Mnumber, OUTCOME2, Population), as.character),
         across(c(B, LL, UL, SE, P, QE, QEp, I2, nstudies, nparticipants, tau2), as.numeric))

meta_results_all <- meta_results_df %>%
  mutate(pFDR = p.adjust(P, method = "fdr")) %>%
  ungroup()

# Load prefix lists from CSVs
M3AgeDisplay <- data.table::fread("M3AgeDisplay.csv")$RISK
M3FIDisplay <- data.table::fread("M3FIDisplay.csv")$RISK

# Function to check if a prefix exists in the list
has_valid_prefix <- function(risk, valid_list) {
  valid_list <- gsub("^g__", "s__", valid_list)  # Replace "g__" with "s__" as we are now looking at species
  valid_list <- gsub("genus", "species", valid_list)  # Replace "genus" with "species" as we are now looking at species
  grepl(paste(valid_list, collapse = "|"), risk) # Check for matches
}

clean_label <- function(varname, table = FALSE) {
  if (is.na(varname)) return(NA_character_)
  
  # Handle species-level variables
  if (stringr::str_starts(varname, "s__")) {
    name <- stringr::str_remove(varname, "^s__")
    toks <- stringr::str_split(name, "_", simplify = TRUE)
    toks <- toks[toks != ""]
    if (length(toks) < 2) {
      return(if (table) name else paste0("italic('", name, "')"))
    }
    
    genus_raw <- toks[1]
    
    if (length(toks) >= 2 && grepl("^\\d+$", toks[2])) {
      genus_key <- paste0(genus_raw, "_", toks[2])   # bv "CAG_217"
      rest <- toks[-c(1,2)]                          # rest bevat nu bv "sp000436335"
    } else {
      genus_key <- genus_raw
      rest <- toks[-1]
    }
    
    # Map placeholder genera to readable genus names with aliases in parentheses
    genus_disp <- dplyr::case_when(
      genus_key %in% c("CAG_217","CAG_273") ~ paste0("Clostridium (sp. CAG:", stringr::str_extract(genus_key, "\\d+"), ")"),
      genus_key == "CAG_83"   ~ "Oscillibacter (sp. CAG:83)",
      genus_key == "ER4"      ~ "Oscillibacter (sp. ER4)",
      genus_key == "SFMI01"   ~ "Christensenellales (sp. SFMI01)",
      genus_key == "UBA1417"  ~ "Acutalibacter (sp. UBA1417)",
      TRUE ~ genus_raw
    )
    
    
    # Identify the species token: first lowercase word (e.g. prausnitzii) or sp + numbers (e.g. sp000765235)
    species_idx <- which(
      grepl("^[a-z][a-z0-9-]*$", rest) |
        grepl("^sp0*\\d+$", rest)
    )
    species_tok <- if (length(species_idx)) rest[species_idx[1]] else if (length(rest)) rest[length(rest)] else NA_character_
    
    # If no valid species token is found, take the last token as a fallback
    if (is.na(species_tok)) species_tok <- rest[length(rest)]
    
    # Keep numeric suffix only if it’s part of a “sp000...” placeholder
    keep_numeric <- grepl("^sp\\d+$", species_tok)
    
    # Construct label
    lab <- if (!is.na(species_tok)) paste(genus_disp, species_tok) else genus_disp
    lab <- gsub("\\s+", " ", lab)
    
    return(if (table) lab else paste0("italic('", lab, "')"))
  }
  
  # Handle genus-level variables
  if (stringr::str_starts(varname, "g__")) {
    name <- stringr::str_remove(varname, "^g__")
    genus <- dplyr::case_when(
      name %in% c("CAG_217", "CAG_273") ~ paste0("Clostridium (sp. CAG:", stringr::str_extract(name, "\\d+"), ")"),
      name == "CAG_83"   ~ "Oscillibacter (sp. CAG:83)",
      name == "ER4"      ~ "Oscillibacter (sp. ER4)",
      name == "SFMI01"   ~ "Christensenellales (sp. SFMI01)",
      name == "UBA1417"  ~ "Acutalibacter (sp. UBA1417)",
      TRUE ~ stringr::str_replace(name, "_.*$", "")
    )
    return(if (table) genus else paste0("italic('", genus, "')"))
  }
  
  # Handle alpha- and beta-diversity labels
  get_type <- function(v) dplyr::case_when(
    stringr::str_detect(v, "genus")   ~ "Genus",
    stringr::str_detect(v, "species") ~ "Species",
    stringr::str_detect(v, "asv")     ~ "ASV",
    stringr::str_detect(v, "feature") ~ "ASV",
    TRUE ~ ""
  )
  
  if (stringr::str_starts(varname, "simpson_e")) {
    type <- get_type(varname)
    return(if (table) paste0("Inverse Simpson (", type, ")") else paste0("Inverse~Simpson~(", type, ")"))
  }
  
  if (stringr::str_starts(varname, "simpson")) {
    type <- get_type(varname)
    return(if (table) paste0("Simpson (", type, ")") else paste0("Simpson~(", type, ")"))
  }
  
  if (stringr::str_starts(varname, "chao1")) {
    type <- get_type(varname)
    return(if (table) paste0("Chao1 (", type, ")") else paste0("Chao1~(", type, ")"))
  }
  
  if (stringr::str_starts(varname, "shannon")) {
    type <- get_type(varname)
    return(if (table) paste0("Shannon (", type, ")") else paste0("Shannon~(", type, ")"))
  }
  
  if (stringr::str_starts(varname, "min_")) {
    beta <- dplyr::case_when(
      stringr::str_detect(varname, "bray") ~ "Bray-Curtis",
      stringr::str_detect(varname, "jacc") ~ "Jaccard",
      stringr::str_detect(varname, "aitch") ~ "Aitchison",
      stringr::str_detect(varname, "uu") ~ "Unweighted UniFrac",
      stringr::str_detect(varname, "wu") ~ "Weighted UniFrac",
      TRUE ~ NA_character_
    )
    type <- get_type(varname)
    return(if (table)
      paste0("Uniqueness based on ", beta, " (", type, ")")
      else
        paste0("Uniqueness~based~on~", beta, "~(", type, ")"))
  }
  
  varname
}





# Generalized function for processing data
process_data <- function(data, outcome, model, valid_list = NULL) {
  data %>%
    mutate(
      diffexpressed = case_when(
        B > 0 & pFDR < 0.05 ~ "UP",
        B < 0 & pFDR < 0.05 ~ "DOWN",
        TRUE ~ "NO"
      ),
      delabel = ifelse(diffexpressed != "NO" & !is.na(RISK),
                       sapply(RISK, clean_label),
                       NA)
      
    ) %>%
    filter(OUTCOME2 == outcome, Mnumber == model)
}

# Process data
metaage <- process_data(meta_results_all, "age", "m3", M3AgeDisplay)
metafi <- process_data(meta_results_all, "fi", "m3", M3FIDisplay)

# List of outcomes to export
outcomes <- c("age", "fi")

# Optional: List of valid prefixes per outcome
valid_prefixes <- list(
  age = M3AgeDisplay,
  fi = M3FIDisplay
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

m3display = full_join(meta_results_age_val, meta_results_fi_val) %>%
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
  )  
# Create separate displays for age and fi outcomes with both original and scaled CIs
m3agedisplayspecies = m3display %>% 
  filter(OUTCOME2 == 'age') %>%
  select(RISK, Studiesdir, nparticipants, B_CI, B_CI_original, pFDR_handig, QE, QEp_handig, I2, tau2, pFDR)

m3fidisplayspecies = m3display %>% 
  filter(OUTCOME2 == 'fi') %>%
  select(RISK, Studiesdir, nparticipants, B_CI, B_CI_original, pFDR_handig, QE, QEp_handig, I2, tau2, pFDR)

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

process_data_all <- function(data, outcome, model, valid_list = NULL) {
  data %>%
    mutate(
      diffexpressed = case_when(
        B > 0 & pFDR < 0.05 & ndirection > 1 ~ "UP",
        B < 0 & pFDR < 0.05 & ndirection > 1 ~ "DOWN",
        TRUE ~ "NO"
      ),
      delabel = ifelse(diffexpressed != "NO" & !is.na(RISK),
                       sapply(RISK, clean_label),
                       NA),
      
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
      
    ) %>%
    filter(OUTCOME2 == outcome, Mnumber == model)
}

# Save output
# Custom filenames per outcome
file_map <- list(
  age = "../Paper/Supplements/SuppTable6.csv",
  fi  = "../Paper/Supplements/SuppTable9.csv"
)

for (outcome in outcomes) {
  data_out <- process_data_all(meta_results_all, outcome, "m3", valid_list = valid_prefixes[[outcome]]) %>%
    filter(Population == "all") %>% 
    mutate(
      Variable = sapply(RISK, clean_label, table = TRUE)
    ) %>%
    select(Outcome, Variable, RISK, Studiesdir, nparticipants, B_CI, pFDR_handig, QE, QEp_handig, I2, tau2) %>%
    rename(., 
           `Microbial Feature` = Variable, 
           `GreenGenes2 name` = RISK,
           `Nstudies (Ndirection)` = Studiesdir,
           Nparticipants = nparticipants,
           `Beta (CI)` = B_CI,
           pFDR = pFDR_handig,
           `QE p-value` = QEp_handig)
  
  if (nrow(data_out) > 0) {
    
    data_out <- data_out %>%
      mutate(across(where(is.list), ~ unlist(.))) %>%
      mutate(across(everything(), as.character))  # <-- key line
    
    filename <- file_map[[outcome]]
    
    write.csv(data_out, filename, row.names = FALSE)
    cat("Saved:", filename, "\n")
    
  } else {
    cat("No valid species found for outcome:", outcome, "\n")
  }
}


# Subset heterogeneous associations
heterogeneous <- full_join(meta_results_age_val, meta_results_fi_val) %>% subset(., QEp < 0.05 | I2 >= 40)

# Unique combinations of RISK and OUTCOME2
heterogeneous_combinations <- unique(heterogeneous[, c("RISK", "OUTCOME2")])

# Define cohorts used in meta-analyses
dataframes_all <- list(FHS = dataframes$FHSshot %>% filter(., Datasplit == 'all'), 
                       MrOS = dataframes$MrOSshot %>% filter(., Datasplit == 'all'), 
                       LL = dataframes$LLshot %>% filter(., Datasplit == 'all'), 
                       SOL = dataframes$SOLshot %>% filter(., Datasplit == 'all')) 

extract_cohort_data <- function(df, cohort_name, risk, outcome) {
  df %>%
    filter(Mnumber == "m3", Variable == risk, Outcome == outcome, Datasplit == 'all') %>%
    mutate(Cohort = paste0(cohort_name, " (N=", N, ")"),
           B = Coefficient) %>%
    select(Cohort, B, LL, UL)
}
# Process each heterogeneous combination
forest_plot_list <- list()

for (i in 1:nrow(heterogeneous_combinations)) {
  risk <- heterogeneous_combinations$RISK[i]
  outcome <- heterogeneous_combinations$OUTCOME2[i]
  outcome2 = ifelse(outcome == 'age', 'Age', ifelse(outcome == 'fi', 'FI', outcome))
  
  # Extract meta-analysis result
  meta_effect <- heterogeneous %>%
    filter(RISK == risk, OUTCOME2 == outcome) %>%
    mutate(Cohort = paste0('Meta-analysis', " (N=", nparticipants, ")")) %>%
    select(Cohort, B, LL, UL)
  
  # Extract data from each cohort
  cohort_results <- do.call(rbind, lapply(names(dataframes_all), function(cohort) {
    extract_cohort_data(dataframes_all[[cohort]], cohort, risk, outcome)
  }))
  
  # Combine meta-analysis with cohort-specific results
  forest_data <- rbind(cohort_results, meta_effect)
  
  # Ensure "Meta-analysis (N=...)" is correctly ordered dynamically
  meta_name <- unique(meta_effect$Cohort)  # Extracts the dynamically generated Meta-analysis name
  
  forest_data$Cohort <- factor(forest_data$Cohort, 
                               levels = c(setdiff(unique(forest_data$Cohort), meta_name), meta_name))
  
  # Create forest plot
  plot <- ggplot(forest_data, aes(y = Cohort, x = B, xmin = LL, xmax = UL)) +
    geom_pointrange(aes(color = Cohort), size = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    theme_minimal(base_size = 15) +
    labs(title = paste("Forest Plot for", risk, "and", outcome2),
         x = "Effect Size (B)", y = "Cohort") +
    theme(legend.position = "none",
          axis.text.y = element_text(face = ifelse(levels(forest_data$Cohort) == meta_name, "bold", "plain")))  # Bold Meta-analysis
  
  
  forest_plot_list[[paste(risk, outcome, sep = "_")]] <- plot
}

# Save plots
for (name in names(forest_plot_list)) {
  ggsave(paste0(name, "_species_forestplot.png"), forest_plot_list[[name]], width = 8, height = 6, units = "in", limitsize = FALSE)
}

# Function for volcano plots
create_volcano_plot <- function(data, title, max_x = NA) {
  filtered_data <- dplyr::filter(data, diffexpressed != "NO")
  
  if (is.na(max_x) && nrow(filtered_data) > 0) {
    extreme_B <- filtered_data$B[which.max(abs(filtered_data$B))]
    max_abs_x <- abs(extreme_B) + 0.1
  } else if (!is.na(max_x)) {
    max_abs_x <- max_x
  } else {
    max_abs_x <- max(abs(data$B), na.rm = TRUE)
  }
  
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


# Create and save plots
plots <- list(
  Age = create_volcano_plot(meta_results_age_m3, "Age", 1), #Filter out s__Lawsonibacter_asaccharolyticus as p < e100 thereby ruining plot
  FI = create_volcano_plot(meta_results_fi_m3, "Frailty Index", 0.05)
)

lapply(names(plots), function(name) {
  ggsave(paste0(name, "species_plot.png"), plots[[name]], width = 15, height = 9, units = "in", limitsize = F)
})

