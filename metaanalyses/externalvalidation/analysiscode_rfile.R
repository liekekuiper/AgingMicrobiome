rm(list = ls())
library(survival)
library(tidyverse)
library(broom)

source("microbiome_utils.R")

# ==================== USER FILLS IN ====================
data_source <- "TSE"          # <-- "TSE"  or  "phyloseq"

# Load the object matching your choice above:
tse_obj    <- readRDS("/path/to/your/tse_object.rds")       # used when data_source == "TSE"
physeq_obj <- readRDS("/path/to/your/phyloseq_object.rds")  # used when data_source == "phyloseq"

metadata    <- "/path/to/your/metadata.txt"
modsfile    <- "/path/to/your/mods_agingmicrobiome.txt"
outsfile    <- "/path/to/your/outs_agingmicrobiome.txt"
cohortname  <- "/path/to/your/cohort.txt"
subsfile    <- "/path/to/your/subsets_agingmicrobiome.txt"
label       <- "16s"           # "16s" or "metagenomics"; controls species-level processing
factors_str <- "NA"            # cohort-specific factors; if >1 use "factor1,factor2"
# =======================================================

# Point to the right object for this run
data_obj <- if (data_source == "TSE") tse_obj else physeq_obj

# Factors
factors <- ifelse(factors_str == "NA", character(), strsplit(factors_str, ",")[[1]])
print(factors)
cat("Imported environments\n")

# Read config files
models   <- readLines(modsfile)
outcomes <- readLines(outsfile)
cohort   <- readLines(cohortname)[1]
subsets  <- readLines(subsfile)
cat("Read files\n")

variables_store <- tolower(colnames(read.delim(metadata)))
cat("Stored variables, will start for loop\n")

# Ensure intermediate output directory exists
dir.create("intermediatefiles", showWarnings = FALSE, recursive = TRUE)

# Run once with strict filter and once with lenient filter
for (lenient in c(FALSE, TRUE)) {
  cat("Starting run with lenient =", lenient, "\n")
  suffix_len <- if (lenient) "_lenient" else ""

  # Fresh results tibble per lenient mode
  results_df <- tibble(
    Datasplit   = character(), Outcome = character(), Variable = character(), Model = character(),
    N           = numeric(),   Ncases  = numeric(),   Coefficient = numeric(), Std.Error = numeric(),
    HR          = numeric(),   LL      = numeric(),   UL          = numeric(),
    t.value     = numeric(),   P       = numeric()
  )

  for (subs in subsets) {
    cat("Performing analyses of", subs, "\n")

    for (out in outcomes) {

      for (model in models) {
        original_model <- model

        if (!(out %in% c("age", "mortality"))) {
          model <- paste0(model, "+age")
        }
        if (grepl("men", subs)) {
          model <- gsub("sex\\+|sex", "", model)
        }
        model_terms <- unlist(strsplit(model, "\\+"))

        cat("Starting to create dataset, for", subs, "-", out, "adjusted for", model, "\n")

        # Pre-process data via unified process() dispatcher
        datafile <- process(
          data_obj      = data_obj,
          metadata_path = metadata,
          model         = model,
          sub           = subs,
          outcome       = out,
          factors       = factors,
          label         = label,
          data_source   = data_source,
          assay_name    = "counts",   # ignored for phyloseq
          lenient       = lenient
        )

        # Sanitise column names
        colnames(datafile) <- gsub("-",          "_", colnames(datafile))
        colnames(datafile) <- gsub(" ",          "_", colnames(datafile))
        colnames(datafile) <- gsub("[^0-9a-zA-Z_]", "", colnames(datafile))

        cat("Created dataset, for", subs, "-", out, "adjusted for", model,
            "— now continuing with analyses\n")

        variables <- setdiff(
          colnames(datafile),
          c(
            variables_store,
            unlist(lapply(model_terms, function(term)
              grep(paste0("^", term, "_"), colnames(datafile), value = TRUE)))
          )
        )

        for (var in variables) {

          if (out != "mortality") {
            # ---- Linear model ----
            formula         <- as.formula(paste(out, "~", model, "+", var))
            datafile[[out]] <- as.numeric(datafile[[out]])
            lm_model        <- lm(formula, data = datafile)
            tidy_res        <- tidy(lm_model)
            row             <- tidy_res[tidy_res$term == var, ]

            if (nrow(row) > 0) {
              results_df <- bind_rows(results_df, tibble(
                Datasplit   = subs,
                Outcome     = out,
                Variable    = var,
                Model       = model,
                N           = nobs(lm_model),
                Ncases      = NA,
                Coefficient = row$estimate,
                Std.Error   = row$std.error,
                HR          = NA,
                LL          = row$estimate - 1.96 * row$std.error,
                UL          = row$estimate + 1.96 * row$std.error,
                t.value     = row$statistic,
                P           = row$p.value
              ))
            }

          } else {
            # ---- Cox PH model ----
            datafile$followup <- datafile$studytime + datafile$age
            covariates        <- unique(trimws(c(var, model_terms)))
            surv_formula      <- as.formula(
              paste("Surv(followup, mortality) ~", paste(covariates, collapse = "+"))
            )

            tryCatch({
              cox_model <- coxph(surv_formula, data = datafile)
              tidy_res  <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
              row       <- tidy_res[tidy_res$term == var, ]

              if (nrow(row) > 0) {
                results_df <- bind_rows(results_df, tibble(
                  Datasplit   = subs,
                  Outcome     = out,
                  Variable    = var,
                  Model       = model,
                  N           = cox_model$n,
                  Ncases      = sum(cox_model$y[, "status"]),
                  Coefficient = row$estimate,
                  Std.Error   = row$std.error,
                  HR          = row$estimate,
                  LL          = row$conf.low,
                  UL          = row$conf.high,
                  t.value     = NA,
                  P           = row$p.value
                ))
              }
            }, error = function(e) {
              cat(sprintf(
                "Failed CoxPH fit for variable '%s' in subset '%s' with outcome '%s': %s\n",
                var, subs, out, e$message
              ))
            })
          }
        } # end variable loop

        # Write intermediate results after each model
        write_csv(
          results_df,
          paste0(
            "intermediatefiles/Results_",
            gsub("[/\\\\]", "_", label), "_",
            subs, "_",
            cohort, "_",
            Sys.Date(),
            suffix_len, ".csv"
          )
        )

        model <- original_model  # restore before next iteration
        cat("Finished analyses of", subs, "for", out, ", adjusted for:", model, "\n")
      } # end model loop

      cat("Finished", out, "analyses for", subs, "\n")
    } # end outcome loop

    cat("Finished analyses of", subs, "\n")
  } # end subset loop

  # Write final results for this lenient mode
  write_csv(
    results_df,
    paste0(
      "Results_",
      cohort, "_",
      gsub("[/\\\\]", "_", label), "_",
      Sys.Date(),
      suffix_len, ".csv"
    )
  )
  cat("Finished run with lenient =", lenient, "\n")
} # end lenient loop
