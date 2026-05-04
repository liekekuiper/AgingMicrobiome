# =============================================================================
# MicrobiomeAging_functions.R
# Functions for microbiome-aging meta-analysis pipeline
# =============================================================================


# -----------------------------------------------------------------------------
# 0. UTILITIES
# -----------------------------------------------------------------------------

`%+%` <- function(a, b) paste0(a, b)

#' Format a p-value for display
#'
#' Values < 0.01 are shown in scientific notation (e.g. 1.23x10^-3),
#' values in [0.045, 0.05) that would round to 0.05 are shown as "<0.05",
#' and all others are shown to two decimal places.
format_value <- function(p) {
  dplyr::case_when(
    p < 0.01                        ~ sub("e", "x10^", sprintf(p, fmt = "%.2e")),
    p < 0.05 & round(p, 2) >= 0.05 ~ "<0.05",
    TRUE                            ~ sprintf(p, fmt = "%.2f")
  )
}

#' Write a supplementary table to the Supplements folder
write_supp_table <- function(df, filename) {
  df %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
    data.table::fwrite(
      file  = file.path("../Paper/Supplements", filename),
      na    = "",
      quote = TRUE
    )
}

#' Test whether a formatted p-value string is significant
is_sig <- function(p) {
  p <- as.character(p)
  grepl("x10|<", p) |
    (!is.na(suppressWarnings(as.numeric(p))) &
       suppressWarnings(as.numeric(p)) < 0.05)
}


# -----------------------------------------------------------------------------
# 1. DATA LOADING & PREPARATION
# -----------------------------------------------------------------------------

#' Expand MrOS dataset to include age_5 and men subgroups
#' (MrOS has no within-study subgroups so these mirror the full dataset)
expand_MrOS <- function(df) {
  dplyr::bind_rows(
    df,
    df %>% dplyr::mutate(Datasplit = "age_5"),
    df %>% dplyr::mutate(Datasplit = "men")
  )
}

#' Load cohort summary files and organise into seq-type sub-lists
#'
#' @param files           Named list of file paths (names encode cohort + seq type).
#' @param add_all         If TRUE, create an "all" entry per cohort.
#' @param expand_fn       Optional function applied to selected datasets.
#' @param expand_datasets Character vector of dataset names to expand.
#' @param mean_cohort_age_list Named list of mean cohort ages (optional).
load_and_prepare_data <- function(files,
                                  add_all              = TRUE,
                                  expand_fn            = NULL,
                                  expand_datasets      = NULL,
                                  mean_cohort_age_list = NULL) {
  
  dataframes <- setNames(lapply(files, data.table::fread, data.table = FALSE), names(files))
  
  if (!is.null(expand_fn) && !is.null(expand_datasets)) {
    for (ds in expand_datasets) {
      if (ds %in% names(dataframes))
        dataframes[[ds]] <- expand_fn(dataframes[[ds]])
    }
  }
  
  dataframes <- mapply(function(df, name) {
    df$cohort <- sub("(16s|shot)$", "", name)
    df
  }, dataframes, names(dataframes), SIMPLIFY = FALSE)
  
  if (!is.null(mean_cohort_age_list)) {
    for (name in names(dataframes))
      dataframes[[name]]$mean_age <- mean_cohort_age_list[[name]]
  }
  
  prefixes  <- unique(sub("(16s|shot)$", "", names(dataframes)))
  data_16s  <- list()
  data_shot <- list()
  data_all  <- list()
  
  for (prefix in prefixes) {
    df_16s  <- dataframes[[paste0(prefix, "16s")]]
    df_shot <- dataframes[[paste0(prefix, "shot")]]
    
    if (!is.null(df_16s))
      data_16s[[paste0(prefix, "16s")]]  <- df_16s  %>% dplyr::mutate(seq_type = "16s")
    
    if (!is.null(df_shot))
      data_shot[[paste0(prefix, "shot")]] <- df_shot %>% dplyr::mutate(seq_type = "shot")
    
    if (add_all) {
      if (!is.null(df_16s) && !is.null(df_shot))
        data_all[[paste0(prefix, "all")]] <- df_shot %>% dplyr::mutate(seq_type = "all")
      else if (!is.null(df_16s))
        data_all[[paste0(prefix, "all")]] <- df_16s  %>% dplyr::mutate(seq_type = "all")
      else if (!is.null(df_shot))
        data_all[[paste0(prefix, "all")]] <- df_shot %>% dplyr::mutate(seq_type = "all")
    }
  }
  
  result        <- c(data_16s = data_16s, data_shot = data_shot, data_all = data_all)
  names(result) <- sub("^[^.]+\\.", "", names(result))
  result
}

#' Derive model identifier (m1–m4) from Model string and drop rows without P
mutate_summ_data <- function(df) {
  df %>%
    dplyr::mutate(
      Mnumber = dplyr::case_when(
        !grepl("statin", Model)                              ~ "m1",
        grepl("statin", Model) & !grepl("bmi", Model)       ~ "m2",
        grepl("bmi", Model)    & !grepl("dietscore", Model) ~ "m3",
        grepl("dietscore", Model)                           ~ "m4",
        TRUE                                                 ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(P))
}


# =============================================================================
# 2. TAXONOMY LABEL CLEANING
# =============================================================================

# Special genus codes → display names (used by both g__ and s__ branches)
resolve_genus <- function(key) {
  dplyr::case_when(
    key %in% c("CAG_217", "CAG_273") ~
      paste0("Clostridium sp. CAG:", stringr::str_extract(key, "\\d+")),
    key == "CAG_83"  ~ "Oscillibacter sp. CAG:83",
    key == "ER4"     ~ "Oscillibacter sp. ER4",
    key == "SFMI01"  ~ "Christensenellales sp. SFMI01",
    key == "UBA1417" ~ "Acutalibacter sp. UBA1417",
    TRUE             ~ stringr::str_replace(key, "_.*$", "")
  )
}

# Render "Genus sp. STRAIN [suffix]": italic genus word, roman remainder.
render_sp_genus <- function(genus_disp, suffix = NULL, table = FALSE) {
  full <- if (!is.null(suffix)) paste(genus_disp, suffix) else genus_disp
  if (table) return(full)
  genus_word <- stringr::str_extract(genus_disp, "^\\S+")
  roman_part <- stringr::str_remove(genus_disp, "^\\S+ ")
  roman_full <- if (!is.null(suffix)) paste(roman_part, suffix) else roman_part
  paste0("italic('", genus_word, "')~'", roman_full, "'")
}

#' Format a microbial variable name as a human-readable label
#'
#' @param varname  Raw variable name (e.g. "g__Faecalibacterium", "s__ER4_sp000765235").
#' @param table    TRUE → plain text; FALSE → plotmath expression string.
clean_label <- function(varname, table = FALSE) {
  if (is.na(varname)) return(NA_character_)
  
  # ── Species: s__Genus_species  or  s__CODE_spXXXXX ──────────────────────────
  if (stringr::str_starts(varname, "s__")) {
    name <- stringr::str_remove(varname, "^s__")
    toks <- stringr::str_split(name, "_", simplify = TRUE)
    toks <- toks[toks != ""]
    
    if (length(toks) < 2)
      return(if (table) name else paste0("italic('", name, "')"))
    
    # Two-token genus key when second token is all digits (e.g. CAG_217)
    if (grepl("^\\d+$", toks[2])) {
      genus_key <- paste0(toks[1], "_", toks[2])
      rest      <- toks[-c(1, 2)]
    } else {
      genus_key <- toks[1]
      rest      <- toks[-1]
    }
    
    genus_disp  <- resolve_genus(genus_key)
    
    # Species token: prefer a lowercase word or spNNNN code
    sp_idx      <- which(grepl("^[a-z][a-z0-9-]*$", rest) | grepl("^sp0*\\d+$", rest))
    species_tok <- if (length(sp_idx)) rest[sp_idx[1]] else rest[length(rest)]
    
    if (table)
      return(gsub("\\s+", " ", paste(genus_disp, species_tok)))
    
    # Placeholder genus ("Genus sp. STRAIN") → italic genus word only, roman rest
    if (grepl(" sp\\. ", genus_disp))
      return(render_sp_genus(genus_disp, species_tok, table = FALSE))
    
    # Numeric strain code (e.g. sp000434635) → italic genus only, roman code
    if (grepl("^sp0*\\d+$", species_tok))
      return(paste0("italic('", genus_disp, "')~'", species_tok, "'"))
    
    # Normal binomial → fully italic
    return(paste0("italic('", genus_disp, " ", species_tok, "')"))
  }
  
  # ── Genus: g__Genus  or  g__CODE ────────────────────────────────────────────
  if (stringr::str_starts(varname, "g__")) {
    name  <- stringr::str_remove(varname, "^g__")
    genus <- resolve_genus(name)
    if (table) return(genus)
    if (grepl(" sp\\. ", genus)) return(render_sp_genus(genus, table = FALSE))
    return(paste0("italic('", genus, "')"))
  }
  
  # ── Alpha-diversity & uniqueness ─────────────────────────────────────────────
  get_level <- function(v) dplyr::case_when(
    stringr::str_detect(v, "genus")       ~ "Genus",
    stringr::str_detect(v, "species")     ~ "Species",
    stringr::str_detect(v, "asv|feature") ~ "ASV",
    TRUE                                  ~ ""
  )
  fmt_div <- function(tbl_name, plt_name, level) {
    if (table) paste0(tbl_name, " (", level, ")")
    else       paste0(plt_name, "~(", level, ")")
  }
  
  if (stringr::str_starts(varname, "simpson_e"))
    return(fmt_div("Inverse Simpson", "Inverse~Simpson", get_level(varname)))
  if (stringr::str_starts(varname, "simpson"))
    return(fmt_div("Simpson", "Simpson", get_level(varname)))
  if (stringr::str_starts(varname, "chao1"))
    return(fmt_div("Chao1", "Chao1", get_level(varname)))
  if (stringr::str_starts(varname, "shannon"))
    return(fmt_div("Shannon", "Shannon", get_level(varname)))
  
  if (stringr::str_starts(varname, "min_")) {
    metric <- dplyr::case_when(
      stringr::str_detect(varname, "bray")  ~ "Bray-Curtis",
      stringr::str_detect(varname, "jacc")  ~ "Jaccard",
      stringr::str_detect(varname, "aitch") ~ "Aitchison",
      stringr::str_detect(varname, "uu")    ~ "Unweighted UniFrac",
      stringr::str_detect(varname, "wu")    ~ "Weighted UniFrac",
      TRUE                                  ~ NA_character_
    )
    level <- get_level(varname)
    return(if (table)
      paste0("Uniqueness based on ", metric, " (", level, ")")
      else
        paste0("Uniqueness~based~on~", metric, "~(", level, ")"))
  }
  
  varname
}

#' Convert a table-format label to a plotmath expression string (for heatmap y-axis)
make_italic_label <- function(lbl) {
  if (grepl("Shannon|Simpson|Chao1|Uniqueness", lbl)) return(lbl)
  
  # "Genus sp. STRAIN ..." → italic genus word, roman remainder
  if (grepl("^\\S+ sp\\. ", lbl)) {
    genus <- regmatches(lbl, regexpr("^\\S+", lbl))
    rest  <- sub("^\\S+ ", "", lbl)
    return(paste0("italic('", genus, "')~'", rest, "'"))
  }
  
  # "Genus spNNNNNN" numeric strain code → italic genus only, roman code
  if (grepl("^\\S+ sp\\d+$", lbl)) {
    genus <- regmatches(lbl, regexpr("^\\S+", lbl))
    code  <- sub("^\\S+ ", "", lbl)
    return(paste0("italic('", genus, "')~'", code, "'"))
  }
  
  # All other taxa (regular binomials) → fully italic
  paste0("italic('", lbl, "')")
}

#' Safely parse a plotmath label into an expression; falls back to plain text
safe_parse <- function(lbl) {
  expr_str <- make_italic_label(lbl)
  if (grepl("italic\\(", expr_str))
    tryCatch(parse(text = expr_str), error = function(e) lbl)
  else
    lbl
}

# -----------------------------------------------------------------------------
# 3. META-ANALYSIS
# -----------------------------------------------------------------------------

#' Check whether a species risk name has a valid genus-level prefix
has_valid_prefix <- function(risk, valid_list) {
  valid_list <- gsub("^g__", "s__", valid_list)
  valid_list <- gsub("genus", "species", valid_list)
  grepl(paste(valid_list, collapse = "|"), risk)
}

#' Run random-effects meta-analyses across all combinations of outcome, model,
#' population, analysis type, and risk factor
#'
#' @param dataframes            Named list of cohort data frames.
#' @param outcome_levels        Character vector of outcomes to analyse.
#' @param mnumber_levels        Character vector of model numbers (e.g. m1–m4).
#' @param population_levels     Character vector of Datasplit values.
#' @param risk_levels_by_outcome Optional named list; if supplied, overrides
#'                              automatic risk discovery per outcome.
#' @param analysis_types        Sequencing types to loop over.
#' @param ncohorts              Maximum expected studies; warns if exceeded.
run_meta_analysis <- function(dataframes,
                              outcome_levels,
                              mnumber_levels,
                              population_levels,
                              risk_levels_by_outcome = NULL,
                              analysis_types         = c("16s", "shot", "all"),
                              ncohorts               = 6) {
  results <- list()
  
  for (outcome in outcome_levels) {
    risk_levels <- if (!is.null(risk_levels_by_outcome))
      risk_levels_by_outcome[[outcome]]
    else
      unique(unlist(lapply(dataframes, function(df) df$Variable)))
    
    for (population in population_levels) {
      for (mnumber in mnumber_levels) {
        for (analysis_type in analysis_types) {
          for (risk in risk_levels) {
            
            all_studies <- lapply(dataframes, function(df) {
              if (risk %in% df$Variable)
                df %>%
                dplyr::filter(seq_type == analysis_type, Variable == risk,
                              Mnumber  == mnumber,       Outcome  == outcome,
                              Datasplit == population) %>%
                dplyr::select(N, Coefficient, Std.Error, P, LL, UL)
            }) %>% dplyr::bind_rows()
            
            if (nrow(all_studies) > ncohorts) {
              message("Too many studies for: ", risk, " | ", outcome, " | ", mnumber)
              next
            }
            if (nrow(all_studies) < 2) next
            
            ndirection <- max(sum(all_studies$Coefficient > 0, na.rm = TRUE),
                              sum(all_studies$Coefficient < 0, na.rm = TRUE))
            
            meta_result <- tryCatch(
              metafor::rma(yi = Coefficient, sei = Std.Error, data = all_studies, method = "REML"),
              error = function(e) {
                message("Falling back to FE: ", risk, " | ", outcome, " | ", mnumber)
                metafor::rma(yi = Coefficient, sei = Std.Error, data = all_studies, method = "FE")
              }
            )
            
            results <- append(results, list(list(
              RISK = risk, Mnumber = mnumber, OUTCOME2 = outcome,
              Population = population, Analysis_Type = analysis_type,
              B = meta_result$b, LL = meta_result$ci.lb, UL = meta_result$ci.ub,
              SE = meta_result$se, P = meta_result$pval,
              QE = meta_result$QE, QEp = meta_result$QEp, I2 = meta_result$I2,
              nstudies = nrow(all_studies), tau2 = meta_result$tau2,
              ndirection = ndirection, nparticipants = sum(all_studies$N)
            )))
          }
        }
      }
    }
  }
  
  do.call(rbind, results) %>%
    as.data.frame() %>%
    dplyr::mutate(
      dplyr::across(c(RISK, Mnumber, OUTCOME2, Population, Analysis_Type), as.character),
      dplyr::across(c(B, LL, UL, SE, P, QE, QEp, I2, nstudies, tau2, ndirection, nparticipants), as.numeric),
      pFDR = p.adjust(P, method = "fdr")
    )
}

#' Run meta-regression with Datasplit as moderator for subgroup analyses
#'
#' @param sig_pairs       Data frame with columns RISK and OUTCOME2 (significant hits).
#' @param dataframes_all  Named list of cohort data frames filtered to Datasplit == "all".
#' @param mnumber         Model to use (default "m3").
#' @param population_groups List of Datasplit group vectors; first element of each
#'                        is used as reference level.
run_meta_regression <- function(sig_pairs,
                                dataframes_all,
                                mnumber          = "m3",
                                population_groups = list(
                                  c("age_3", "age_2", "age_1", "age_4", "age_5"),
                                  c("women", "men")
                                )) {
  subgroup_results <- list()
  
  for (i in seq_len(nrow(sig_pairs))) {
    risk    <- sig_pairs$RISK[i]
    outcome <- sig_pairs$OUTCOME2[i]
    
    for (pop_group in population_groups) {
      all_studies <- lapply(dataframes_all, function(df) {
        df_f <- df %>%
          dplyr::filter(Mnumber == mnumber, Outcome == outcome, Datasplit %in% pop_group)
        if (risk %in% df_f$Variable)
          df_f %>% dplyr::filter(Variable == risk) %>%
          dplyr::select(N, Coefficient, Std.Error, P, LL, UL, Datasplit, cohort)
        else NULL
      }) %>% dplyr::bind_rows()
      
      if (nrow(all_studies) < 2) next
      
      ref <- pop_group[1]
      all_studies <- all_studies %>%
        dplyr::mutate(Datasplit = factor(Datasplit, levels = c(ref, setdiff(pop_group, ref))))
      
      meta_result <- tryCatch(
        metafor::rma(yi = Coefficient, sei = Std.Error, mods = ~ Datasplit,
                     data = all_studies, method = "REML"),
        error = function(e) {
          message("Falling back to FE: ", risk, " | ", outcome)
          metafor::rma(yi = Coefficient, sei = Std.Error, mods = ~ Datasplit,
                       data = all_studies, method = "FE")
        }
      )
      
      estimates           <- as.data.frame(meta_result$b)
      colnames(estimates) <- "B"
      estimates$SE        <- meta_result$se
      estimates$P         <- meta_result$pval
      estimates$LL        <- meta_result$ci.lb
      estimates$UL        <- meta_result$ci.ub
      estimates$Datasplit <- rownames(estimates)
      
      for (j in seq_len(nrow(estimates))) {
        subgroup_name <- gsub("Datasplit", "", estimates$Datasplit[j])
        
        if (subgroup_name == "intrcpt") {
          subgroup_name          <- ref
          nstudies_sub           <- length(unique(all_studies$cohort))
          nparticipants_sub      <- sum(all_studies$N, na.rm = TRUE)
          ndirection_sub         <- "-"
        } else {
          sub_data               <- dplyr::filter(all_studies, Datasplit == subgroup_name)
          nstudies_sub           <- nrow(sub_data)
          nparticipants_sub      <- sum(sub_data$N, na.rm = TRUE)
          ndirection_sub         <- max(sum(sub_data$Coefficient > 0, na.rm = TRUE),
                                        sum(sub_data$Coefficient < 0, na.rm = TRUE))
        }
        
        subgroup_results <- append(subgroup_results, list(list(
          RISK             = risk,
          OUTCOME2         = outcome,
          Mnumber          = mnumber,
          Population_Group = subgroup_name,
          B                = estimates$B[j],
          SE               = estimates$SE[j],
          P                = estimates$P[j],
          LL               = estimates$LL[j],
          UL               = estimates$UL[j],
          QE               = meta_result$QE,
          QEp              = meta_result$QEp,
          I2               = meta_result$I2,
          nstudies         = nstudies_sub,
          tau2             = meta_result$tau2,
          nparticipants    = nparticipants_sub,
          QM               = meta_result$QM,
          QMp              = meta_result$QMp,
          ndirection       = ndirection_sub,
          R2               = meta_result$R2
        )))
      }
    }
  }
  
  meta_regression_df <- do.call(rbind, subgroup_results) %>%
    as.data.frame() %>%
    dplyr::mutate(P = as.numeric(P), QMp = as.numeric(QMp))
  
  # FDR-correct at unique-test level then join back
  unique_tests <- meta_regression_df %>%
    dplyr::select(RISK, OUTCOME2, Mnumber, I2, tau2, QM, QE, QMp, QEp) %>%
    dplyr::distinct() %>%
    dplyr::mutate(QMpFDR = p.adjust(QMp, method = "fdr"),
                  QEpFDR = p.adjust(QEp, method = "fdr"))
  
  dplyr::left_join(meta_regression_df, unique_tests) %>%
    dplyr::mutate(pFDR = p.adjust(P, method = "fdr")) %>%
    dplyr::group_by(RISK, OUTCOME2, Mnumber, I2, tau2) %>%
    dplyr::mutate(
      show_tests     = dplyr::row_number() == 2,
      QM_display     = ifelse(show_tests, QM,     NA),
      QMpFDR_display = ifelse(show_tests, QMpFDR, NA),
      QE_display     = ifelse(show_tests, QE,     NA),
      QEpFDR_display = ifelse(show_tests, QEpFDR, NA),
      R2_display     = dplyr::case_when(
        QMpFDR < 0.05  & show_tests ~ sprintf("%.2f", as.numeric(R2)),
        QMpFDR >= 0.05 & show_tests ~ "-",
        !show_tests                  ~ NA_character_
      )
    ) %>%
    dplyr::ungroup()
}


# -----------------------------------------------------------------------------
# 4. DATA PROCESSING FOR FIGURES & TABLES
# -----------------------------------------------------------------------------

#' Extract cohort-level estimates for a single RISK × outcome combination
extract_cohort_data <- function(df, cohort_name, risk, outcome) {
  df_sub <- df %>%
    dplyr::filter(Mnumber == "m3", Variable == risk, Outcome == outcome, Datasplit == "all")
  if (nrow(df_sub) == 0) return(NULL)
  df_sub %>%
    dplyr::mutate(Cohort = paste0(cohort_name, " (N=", N, ")"), B = Coefficient) %>%
    dplyr::select(Cohort, B, LL, UL)
}

#' Classify associations as UP/DOWN and generate volcano plot labels
process_data <- function(data, outcome, model, valid_list = NULL) {
  data %>%
    dplyr::mutate(
      diffexpressed = dplyr::case_when(
        B > 0 & pFDR < 0.05 ~ "UP",
        B < 0 & pFDR < 0.05 ~ "DOWN",
        TRUE                 ~ "NO"
      ),
      delabel = ifelse(diffexpressed != "NO" & !is.na(RISK),
                       sapply(RISK, clean_label), NA)
    ) %>%
    dplyr::filter(OUTCOME2 == outcome, Mnumber == model)
}

#' Add all display columns needed for supplementary tables
#'
#' Optionally filters to a specific outcome and/or model number.
process_data_suppl <- function(data, outcome = NULL, model = NULL) {
  data %>%
    dplyr::mutate(
      diffexpressed = dplyr::case_when(
        B > 0 & pFDR < 0.05 & ndirection > 1 ~ "UP",
        B < 0 & pFDR < 0.05 & ndirection > 1 ~ "DOWN",
        TRUE                                   ~ "NO"
      ),
      delabel      = ifelse(diffexpressed != "NO" & !is.na(RISK),
                            sapply(RISK, clean_label), NA),
      Variable     = sapply(RISK, clean_label, table = TRUE),
      B_CI         = dplyr::case_when(
        OUTCOME2 == "fi" ~
          paste0(sprintf(B,  fmt = "%.4f"), " (", sprintf(LL, fmt = "%.4f"), "; ", sprintf(UL, fmt = "%.4f"), ")"),
        TRUE ~
          paste0(sprintf(B,  fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), "; ", sprintf(UL, fmt = "%.2f"), ")")
      ),
      pFDR_display = format_value(pFDR),
      QE           = sprintf(QE,   fmt = "%.2f"),
      QEp_display  = format_value(QEp),
      I2           = sprintf(I2,   fmt = "%.2f"),
      tau2         = sprintf(tau2, fmt = "%.2f"),
      Outcome      = dplyr::case_when(OUTCOME2 == "age" ~ "Chronological age",
                                      OUTCOME2 == "fi"  ~ "Frailty index"),
      Studiesdir   = paste0(nstudies, " (", ndirection, ")")
    ) %>%
    { if (!is.null(outcome)) dplyr::filter(., OUTCOME2 == outcome) else . } %>%
    { if (!is.null(model))   dplyr::filter(., Mnumber  == model)   else . }
}

#' Standard select + rename for all supplementary tables
prepare_supp_table <- function(data) {
  data %>%
    dplyr::select(Outcome, Variable, RISK, Studiesdir, nparticipants,
                  B_CI, pFDR_display, QE, QEp_display, I2, tau2) %>%
    dplyr::rename(
      `Microbial Feature`     = Variable,
      `GreenGenes2 name`      = RISK,
      `Nstudies (Ndirection)` = Studiesdir,
      Nparticipants           = nparticipants,
      `Beta (CI)`             = B_CI,
      pFDR                    = pFDR_display,
      `QE p-value`            = QEp_display
    )
}

#' Build the display version of significant meta-analysis results
#'
#' Scales B/LL/UL, formats CIs and p-values, adds Variable label.
#' Pass a custom \code{beta_scale_fn(x, risk, outcome)} for species-level scaling.
build_m3display <- function(data, fi_scale = 100, beta_scale_fn = NULL) {
  data %>%
    dplyr::mutate(
      B_original  = B,  LL_original = LL,  UL_original = UL,
      dplyr::across(
        c(B, LL, UL),
        ~ if (!is.null(beta_scale_fn)) beta_scale_fn(.x, RISK, OUTCOME2)
        else dplyr::case_when(OUTCOME2 == "fi" ~ .x * fi_scale, TRUE ~ .x)
      ),
      Variable      = sapply(RISK, clean_label, table = TRUE),
      B_CI          = paste0(sprintf(B,  fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), "; ", sprintf(UL, fmt = "%.2f"), ")"),
      B_CI_original = dplyr::case_when(
        OUTCOME2 == "fi" ~
          paste0(sprintf(B_original,  fmt = "%.4f"), " (", sprintf(LL_original, fmt = "%.4f"), "; ", sprintf(UL_original, fmt = "%.4f"), ")"),
        TRUE ~
          paste0(sprintf(B_original,  fmt = "%.2f"), " (", sprintf(LL_original, fmt = "%.2f"), "; ", sprintf(UL_original, fmt = "%.2f"), ")")
      ),
      pFDR_display = format_value(pFDR),
      QE           = sprintf(QE,   fmt = "%.2f"),
      QEp_display  = format_value(QEp),
      I2           = sprintf(I2,   fmt = "%.2f"),
      tau2         = sprintf(tau2, fmt = "%.2f"),
      Studiesdir   = paste0(nstudies, " (", ndirection, ")")
    )
}

#' Build a discovery + validation merged table for one outcome
#'
#' @param charge_df  Discovery data frame (from M3*Display.csv).
#' @param val_mg     Validation data (metagenomics cohort).
#' @param val_16s    Validation data (16S cohort).
#' @param charge_out Outcome string used to filter charge_df (unused here; kept for API clarity).
#' @param val_out    Outcome string used to filter val_mg / val_16s.
#' @param level      "Genus" or "Species".
build_table <- function(charge_df, val_mg, val_16s, charge_out, val_out, level) {
  val_mg_sub  <- if ("Outcome" %in% names(val_mg))  val_mg[Outcome  == val_out] else val_mg[0]
  val_16s_sub <- if ("Outcome" %in% names(val_16s)) val_16s[Outcome == val_out] else val_16s[0]
  
  val_cols <- c("Variable", "N", "Ncases", "B_CI", "pFDR_display", "abprev")
  if (val_out == "age") val_cols <- setdiff(val_cols, "Ncases")
  
  prep_subset <- function(sub_df, cols, suffix) {
    if (nrow(sub_df) > 0) {
      dt          <- sub_df[, ..cols]
      rename_cols <- setdiff(cols, "Variable")
      data.table::setnames(dt, rename_cols, paste0(rename_cols, "_", suffix))
    } else {
      dt <- data.table::as.data.table(matrix(NA_character_, nrow = 0, ncol = length(cols)))
      data.table::setnames(dt, cols)
      rename_cols <- setdiff(cols, "Variable")
      data.table::setnames(dt, rename_cols, paste0(rename_cols, "_", suffix))
      dt[, Variable := as.character(Variable)]
    }
    dt
  }
  
  vm  <- prep_subset(val_mg_sub,  val_cols, "MG")
  v16 <- prep_subset(val_16s_sub, val_cols, "16S")
  
  res <- merge(charge_df, vm,  by.x = "RISK", by.y = "Variable", all.x = TRUE)
  res <- merge(res,       v16, by.x = "RISK", by.y = "Variable", all.x = TRUE)
  
  res[, Variable := sapply(RISK, clean_label, table = TRUE)]
  res[, Level    := level]
  
  na_label <- "Not present at prevalence and abundance thresholds"
  
  if ("abprev_MG" %in% names(res)) {
    res[, abprev_MG := as.character(abprev_MG)]
    res[is.na(abprev_MG), abprev_MG := na_label]
  }
  
  if ("abprev_16S" %in% names(res)) {
    res[, abprev_16S := as.character(abprev_16S)]
    if (level == "Genus") {
      res[is.na(abprev_16S), abprev_16S := na_label]
    } else {
      cols_16s <- grep("_16S$", names(res), value = TRUE)
      res[, (cols_16s) := NA_character_]
    }
  }
  
  charge_cols  <- c("Variable", "Level", "Studiesdir", "nparticipants",
                    "B_CI", "pFDR_display", "QE", "QEp_display", "I2")
  val_order_mg  <- grep("_MG$",  names(res), value = TRUE)
  val_order_16s <- grep("_16S$", names(res), value = TRUE)
  
  res[, c(charge_cols, val_order_mg, val_order_16s, "pFDR"), with = FALSE]
}

#' Sort a merged discovery + validation table by validation status and pFDR
apply_custom_sort <- function(dt) {
  dt[, SortKey := 99L]
  dt[Level == "Genus"   & abprev_MG == "strict" & abprev_16S == "strict", SortKey := 1L]
  dt[Level == "Genus"   & abprev_MG == "strict"  & SortKey == 99,         SortKey := 2L]
  dt[Level == "Genus"   & abprev_16S == "strict" & SortKey == 99,         SortKey := 3L]
  dt[Level == "Genus"   & abprev_MG == "lenient" & SortKey == 99,         SortKey := 4L]
  dt[Level == "Genus"   & abprev_16S == "lenient"& SortKey == 99,         SortKey := 5L]
  dt[Level == "Genus"   & SortKey == 99,                                  SortKey := 6L]
  dt[Level == "Species" & abprev_MG == "strict",                          SortKey := 7L]
  dt[Level == "Species" & abprev_MG == "lenient",                         SortKey := 8L]
  dt[Level == "Species" & SortKey == 99,                                  SortKey := 9L]
  
  data.table::setorder(dt, SortKey, pFDR)
  dt[, c("SortKey", "Level", "pFDR") := NULL]
  dt
}

#' Set Feature as a reversed factor (bottom-to-top ordering on y-axis)
set_feat_levels <- function(df) {
  df[, Feature := factor(Feature, levels = rev(unique(as.character(Feature))))]
  df
}

#' Build a binary presence/absence data frame for upset plots
get_risk <- function(df, risk_level) {
  df <- df %>% dplyr::select(Variable) %>% dplyr::distinct() %>%
    dplyr::filter(grepl("^g__", Variable))
  df[[risk_level]] <- 1L
  df
}


# -----------------------------------------------------------------------------
# 5. PLOTTING
# -----------------------------------------------------------------------------

#' Forest plot coloured by sequencing type across significant associations
#'
#' @param meta_results  Full meta-analysis results data frame.
#' @param m3sig         Significant hits (used for semi-join filter).
#' @param outcome       Outcome string ("age" or "fi").
#' @param display_name  Human-readable outcome name for plot title.
#' @param filter        If TRUE, restrict to heterogeneous associations only.
#' @param heterogeneous Data frame of heterogeneous associations (required if filter = TRUE).
create_forest_plot <- function(meta_results, m3sig, outcome, display_name,
                               filter = FALSE, heterogeneous = NULL) {
  filtered_meta <- meta_results %>%
    dplyr::filter(OUTCOME2 == outcome, Population == "all", Mnumber == "m3") %>%
    dplyr::semi_join(m3sig[, c("RISK", "OUTCOME2")], by = c("RISK", "OUTCOME2")) %>%
    dplyr::mutate(
      dplyr::across(c(B, LL, UL),
                    ~ ifelse(grepl("min_|shannon|simpson|chao1", RISK), .x * 0.1, .x)),
      Analysis_Type = factor(Analysis_Type, levels = c("shot", "16s", "all"))
    )
  
  if (filter && !is.null(heterogeneous))
    filtered_meta <- filtered_meta %>%
      dplyr::filter(RISK %in% dplyr::filter(heterogeneous, OUTCOME2 == outcome)$RISK)
  
  filtered_meta <- filtered_meta %>%
    dplyr::arrange(RISK, Analysis_Type) %>%
    dplyr::mutate(
      y_position = dplyr::row_number(),
      y_label    = ifelse(Analysis_Type == "all",
                          sapply(RISK, clean_label, table = FALSE), NA_character_)
    )
  
  ggplot2::ggplot(filtered_meta,
                  ggplot2::aes(x = B, y = factor(y_position), xmin = LL, xmax = UL, color = Analysis_Type)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbarh(height = 0.2) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    ggplot2::labs(x = "Effect Estimate", y = NULL,
                  title = paste("Forest Plot for", display_name)) +
    ggplot2::scale_y_discrete(
      labels = function(x)
        parse(text = filtered_meta$y_label[match(x, as.character(filtered_meta$y_position))])
    ) +
    ggplot2::scale_color_manual(
      name   = "Included datasets",
      values = c(shot = "#648fff", `16s` = "#ffb000", all = "#dc267f"),
      labels = c(shot = "Metagenomics", `16s` = "16S",
                 all  = "All (including metagenomics\n if available)"),
      breaks = rev(levels(filtered_meta$Analysis_Type))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y      = ggplot2::element_text(size = 12),
      axis.ticks.y     = ggplot2::element_blank(),
      legend.position  = "bottom",
      text             = ggplot2::element_text(size = 16)
    )
}

#' Per-cohort forest plots for heterogeneous associations, with optional patchwork export
#'
#' @param meta_results    Data frame already filtered to significant hits (will be
#'                        internally filtered to QEp < 0.05 | I2 >= 40).
#' @param dataframes_all  Named list of cohort data frames (Datasplit == "all").
#' @param cohort_colours  Optional named colour vector for CohortShort values.
#' @param save_individual Save one PNG per association.
#' @param save_patchwork  Save combined PNG + PDF per outcome.
#' @param save_dir        Output directory.
create_heterogeneity_forest_plots <- function(meta_results,
                                              dataframes_all,
                                              cohort_colours  = NULL,
                                              save_individual = TRUE,
                                              save_patchwork  = TRUE,
                                              save_dir        = ".") {
  heterogeneous             <- meta_results %>% dplyr::filter(QEp < 0.05 | I2 >= 40)
  heterogeneous_combinations <- unique(heterogeneous[, c("RISK", "OUTCOME2")])
  
  forest_plot_list <- list()
  
  for (i in seq_len(nrow(heterogeneous_combinations))) {
    risk    <- heterogeneous_combinations$RISK[i]
    outcome <- heterogeneous_combinations$OUTCOME2[i]
    
    meta_effect <- heterogeneous %>%
      dplyr::filter(RISK == risk, OUTCOME2 == outcome) %>%
      dplyr::mutate(Cohort = paste0("Meta-analysis (N=", nparticipants, ")")) %>%
      dplyr::select(Cohort, B, LL, UL)
    
    cohort_results <- do.call(rbind, Filter(Negate(is.null),
                                            lapply(names(dataframes_all), function(cohort)
                                              extract_cohort_data(dataframes_all[[cohort]], cohort, risk, outcome))
    ))
    
    if (is.null(cohort_results) || nrow(cohort_results) == 0) next
    
    meta_name  <- unique(meta_effect$Cohort)
    forest_data <- rbind(cohort_results, meta_effect) %>%
      dplyr::mutate(
        Cohort      = factor(Cohort, levels = c(setdiff(unique(Cohort), meta_name), meta_name)),
        CohortShort = ifelse(grepl("^Meta-analysis", Cohort), "Meta-analysis",
                             sub(" \\(N=.*\\)$", "", as.character(Cohort)))
      )
    
    label_expr <- tryCatch(parse(text = clean_label(risk, table = FALSE))[[1]],
                           error = function(e) risk)
    
    p <- ggplot2::ggplot(forest_data,
                         ggplot2::aes(y = Cohort, x = B, xmin = LL, xmax = UL)) +
      ggplot2::geom_pointrange(ggplot2::aes(color = CohortShort), size = 1.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::labs(title = label_expr, x = "Effect Size (B)", y = "Cohort") +
      ggplot2::theme(
        legend.position = "none",
        axis.text.y     = ggplot2::element_text(
          face = ifelse(levels(forest_data$Cohort) == meta_name, "bold", "plain"))
      )
    
    if (!is.null(cohort_colours))
      p <- p + ggplot2::scale_color_manual(values = cohort_colours)
    
    key <- paste(risk, outcome, sep = "_")
    forest_plot_list[[key]] <- p
    
    if (save_individual)
      ggplot2::ggsave(file.path(save_dir, paste0(key, "_forestplot.png")),
                      p, width = 8, height = 6, units = "in", limitsize = FALSE)
  }
  
  if (save_patchwork) {
    plots_by_outcome <- split(forest_plot_list, sub(".*_", "", names(forest_plot_list)))
    for (outcome in names(plots_by_outcome)) {
      plots <- plots_by_outcome[[outcome]]
      if (length(plots) == 0) next
      combined <- patchwork::wrap_plots(plots, ncol = 2, tag_level = "new") +
        patchwork::plot_annotation(tag_levels = "a")
      h <- 5 * ceiling(length(plots) / 2)
      for (ext in c("png", "pdf"))
        ggplot2::ggsave(file.path(save_dir, paste0("patchwork_", outcome, ".", ext)),
                        combined, width = 15, height = h, units = "in")
    }
  }
  
  invisible(forest_plot_list)
}

#' Scatter plot of cohort-level coefficients across age/sex subgroups
#'
#' @param dataframes_all    Named list of cohort data frames.
#' @param variable          RISK variable name.
#' @param outcome           Outcome string.
#' @param datasplit_groups  Datasplit values to include.
#' @param filename          Output file path.
#' @param line              If TRUE, overlay a smoothing line.
plot_microbiome_effect <- function(dataframes_all, variable, outcome,
                                   datasplit_groups, filename, line = TRUE) {
  df_results <- dplyr::bind_rows(dataframes_all) %>%
    dplyr::filter(Mnumber == "m3", Outcome == outcome, Datasplit %in% datasplit_groups)
  
  if (!variable %in% df_results$Variable) return(NULL)
  
  df_results <- df_results %>%
    dplyr::filter(Variable == variable) %>%
    dplyr::mutate(
      age_use = dplyr::case_when(
        !grepl("age_", Datasplit) ~ mean_age,
        Datasplit == "age_1"      ~ 30,
        Datasplit == "age_2"      ~ 45,
        Datasplit == "age_3"      ~ 55,
        Datasplit == "age_4"      ~ 65,
        Datasplit == "age_5"      ~ 75
      )
    )
  
  title_label <- parse(text = clean_label(variable, table = FALSE))
  
  baseplot <- ggplot2::ggplot(df_results,
                              ggplot2::aes(x = age_use, y = Coefficient, size = N, color = cohort)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_x_continuous(
      breaks = c(30, 45, 55, 65, 75),
      labels = c("<40", "40\u201350", "50\u201360", "60\u201370", "70<")
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60") +
    ggplot2::xlab("Age category") +
    ggplot2::ylab(paste(stringr::str_to_title(outcome), "coefficient")) +
    ggplot2::labs(colour = "Cohort", size = "N participants") +
    ggplot2::ggtitle(label = title_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  
  if (line) {
    plot <- if (grepl("age_", datasplit_groups[1]))
      baseplot + ggplot2::stat_smooth(method = "lm", formula = y ~ poly(x, 3),
                                      se = TRUE, linewidth = 1, color = "black")
    else
      baseplot + ggplot2::stat_smooth(
        ggplot2::aes(group = Datasplit, color = Datasplit),
        method = "lm", formula = y ~ x, se = TRUE, linewidth = 1)
  } else {
    plot <- baseplot
  }
  
  ggplot2::ggsave(filename, plot, dpi = 1200, width = 7, height = 5,
                  units = "in", limitsize = FALSE)
  list(plot = plot, results = df_results)
}


# -----------------------------------------------------------------------------
# 6. HEATMAP HELPERS
# -----------------------------------------------------------------------------

#' Diverging fill scale centred on zero (used in heatmap panels)
#'
#' @param fill_abs Absolute maximum for scale limits.
#' @param legend   Whether to show the colour bar legend.
make_fill <- function(fill_abs, legend = FALSE) {
  ggplot2::scale_fill_gradient2(
    low      = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-fill_abs, fill_abs),
    na.value = "transparent", name = "Effect\nsize",
    guide    = if (legend)
      ggplot2::guide_colorbar(
        barwidth    = ggplot2::unit(0.45, "cm"),
        barheight   = ggplot2::unit(4,    "cm"),
        ticks.colour = "grey40",
        frame.colour = "grey40"
      )
    else "none"
  )
}

#' Build one heatmap panel (alpha, uniqueness, genus or species)
#'
#' @param sub_dt    data.table with columns Feature, Outcome, Effect, Annotation, FeatureGroup.
#' @param strip_col Hex colour for the facet strip background.
#' @param fill_abs  Passed to make_fill().
#' @param show_x    Show x-axis labels (set TRUE for bottom panels).
#' @param show_leg  Show colour bar legend (set TRUE for the last panel).
make_panel <- function(sub_dt, strip_col, fill_abs, show_x = FALSE, show_leg = FALSE) {
  ggplot2::ggplot(sub_dt, ggplot2::aes(x = Outcome, y = Feature, fill = Effect)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = Annotation),
                       size = 4.2, vjust = 0.55, colour = "grey15") +
    ggplot2::facet_wrap(~ FeatureGroup, strip.position = "top") +
    ggplot2::scale_y_discrete(labels = function(x) sapply(x, safe_parse)) +
    make_fill(fill_abs, legend = show_leg) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(fill = strip_col, colour = "grey75", linewidth = 0.4),
      strip.text       = ggplot2::element_text(face = "bold", size = 9, colour = "grey20",
                                               margin = ggplot2::margin(t = 3, b = 3)),
      axis.text.y      = ggplot2::element_text(size = 7.5, colour = "grey25"),
      axis.text.x      = if (show_x)
        ggplot2::element_text(size = 9, colour = "grey20", margin = ggplot2::margin(t = 3))
      else ggplot2::element_blank(),
      axis.ticks.x     = if (show_x)
        ggplot2::element_line(colour = "grey70", linewidth = 0.3)
      else ggplot2::element_blank(),
      axis.ticks.y      = ggplot2::element_line(colour = "grey70", linewidth = 0.3),
      axis.ticks.length = ggplot2::unit(2, "pt"),
      legend.position   = if (show_leg) "right" else "none",
      legend.title      = ggplot2::element_text(size = 8,   colour = "grey30"),
      legend.text       = ggplot2::element_text(size = 7.5),
      plot.margin       = ggplot2::margin(t = 2, r = 4, b = 2, l = 4)
    )
}
