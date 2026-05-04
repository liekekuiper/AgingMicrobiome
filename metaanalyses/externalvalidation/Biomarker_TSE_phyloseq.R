rm(list = ls())

# ============================================================
#  USER INPUT
# ============================================================
data_source  <- "TSE"          # <-- "TSE"  or  "phyloseq"

# Paths / identifiers
tse_path     <- "GG2_TSE.rds"          # used when data_source == "TSE"
ps_path      <- "GG2_physeq.rds"       # used when data_source == "phyloseq"
metafilepath <- "metadata_agingmicrobiome.txt"
cohort       <- "FINRISK"
type         <- "16S"                  # e.g. "16S" or "WGS"
betas_csv    <- "~/Microbiome/M3_B_FIDisplay.csv"
# ============================================================

cat(sprintf("Starting Biomarker script  |  data_source = %s\n", data_source))

# ---- Load packages always needed for the microbiome step ----
if (!require("vegan", quietly = TRUE)) install.packages("vegan")

# ============================================================
#  DATA LOADING  (TSE  vs  phyloseq)
# ============================================================

if (data_source == "TSE") {

  # ---- TreeSummarizedExperiment ----
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    BiocManager::install("SummarizedExperiment")
  library(SummarizedExperiment)

  tse <- readRDS(tse_path)

  # Counts (features x samples)
  if ("counts" %in% assayNames(tse)) {
    counts <- assays(tse)[["counts"]]
  } else {
    counts <- assay(tse, 1)
  }

  # Taxonomy table
  tx  <- as.data.frame(rowData(tse), stringsAsFactors = FALSE)
  fid <- rownames(tse)
  if (is.null(fid) || anyNA(fid) || any(fid == "")) {
    cand_id_cols <- c("FeatureID","feature_id","ASV","asv","OTU","otu","ID","id")
    hit <- cand_id_cols[cand_id_cols %in% colnames(tx)][1]
    if (!is.na(hit)) fid <- as.character(tx[[hit]]) else
      stop("Could not determine FeatureID.")
  }
  rownames(counts) <- fid

  # Helper: pick the first matching column or return empty strings
  pick_rank <- function(df, choices, n) {
    col <- choices[choices %in% colnames(df)][1]
    if (!is.na(col)) {
      v <- as.character(df[[col]]); v[is.na(v)] <- ""
      v <- sub("^[a-z]__", "", v)
    } else v <- rep("", n)
    v
  }
  pref    <- function(x, pfx) paste0(pfx, x)
  n_feat  <- nrow(tx)

  taxon_string <- paste(
    pref(pick_rank(tx, c("Kingdom","kingdom","Domain","domain","Superkingdom","superkingdom","Rank1"), n_feat), "k__"),
    pref(pick_rank(tx, c("Phylum","phylum","Rank2"),  n_feat), "p__"),
    pref(pick_rank(tx, c("Class","class","Rank3"),    n_feat), "c__"),
    pref(pick_rank(tx, c("Order","order","Rank4"),    n_feat), "o__"),
    pref(pick_rank(tx, c("Family","family","Rank5"),  n_feat), "f__"),
    pref(pick_rank(tx, c("Genus","genus","Rank6"),    n_feat), "g__"),
    pref(pick_rank(tx, c("Species","species","Rank7"),n_feat), "s__"),
    sep = "; "
  )

  taxonomy_table <- data.frame(FeatureID = fid, Taxon = taxon_string,
                                stringsAsFactors = FALSE, row.names = NULL)
  stopifnot(nrow(counts) == nrow(taxonomy_table),
            identical(rownames(counts), taxonomy_table$FeatureID))

} else if (data_source == "phyloseq") {

  # ---- phyloseq ----
  if (!requireNamespace("phyloseq", quietly = TRUE)) install.packages("phyloseq")
  library(phyloseq)

  ps <- readRDS(ps_path)
  print(ps)

  tx <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

  pref <- function(x, pfx) {
    x <- as.character(x); x[is.na(x)] <- ""; paste0(pfx, x)
  }

  taxon_string <- paste(
    pref(tx[, intersect("Kingdom", colnames(tx))], "k__"),
    pref(tx[, intersect("Phylum",  colnames(tx))], "p__"),
    pref(tx[, intersect("Class",   colnames(tx))], "c__"),
    pref(tx[, intersect("Order",   colnames(tx))], "o__"),
    pref(tx[, intersect("Family",  colnames(tx))], "f__"),
    pref(tx[, intersect("Genus",   colnames(tx))], "g__"),
    pref(tx[, intersect("Species", colnames(tx))], "s__"),
    sep = "; "
  )
  taxonomy_table <- data.frame(FeatureID = rownames(tx), Taxon = taxon_string,
                                stringsAsFactors = FALSE)

  counts <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) counts <- t(counts)

} else {
  stop(sprintf("Unknown data_source '%s'. Set to \"TSE\" or \"phyloseq\".", data_source))
}

cat("Number of unique taxa:", nrow(taxonomy_table), "\n")

# ============================================================
#  SHARED PROCESSING
# ============================================================

# Metadata
metadata <- read.table(metafilepath, header = TRUE, sep = "\t")
names(metadata) <- tolower(names(metadata))
metadata <- metadata[!is.na(metadata$mortality) & !is.na(metadata$sampleid), ]
cat("Metadata rows after NA filter:", nrow(metadata), "\n")

# Genus key
taxonomy_table$genus <- sub("; s__.*", "", taxonomy_table[, 2])

# Reorder counts to taxonomy order
counts <- counts[taxonomy_table[, 1], , drop = FALSE]

# Collapse to genera
counts_genera <- aggregate(counts ~ factor(taxonomy_table$genus), FUN = sum)
rownames(counts_genera) <- counts_genera[, 1]
counts_genera <- counts_genera[, -1]
counts_genera <- t(counts_genera)

sample_ids    <- as.character(metadata$sampleid)
overlap       <- intersect(sample_ids, rownames(counts_genera))
counts_genera <- counts_genera[overlap, , drop = FALSE]
cat("Overlapping IDs:", length(overlap), "\n")
cat("Dimensions counts_genera:", dim(counts_genera), "\n")

# Abundance & presence filtering
abundance_genera      <- sweep(counts_genera, 1, rowSums(counts_genera), "/")
filter                <- colSums(abundance_genera > 0) > 0.1 * nrow(abundance_genera) &
                         apply(abundance_genera, 2, \(x) mean(x[x > 0]) >= 0.0001)
counts_genera_filtered <- counts_genera[, filter]

genera_clr <- vegan::decostand(counts_genera_filtered, method = "clr", pseudocount = 1)

# Betas & risk score
betas         <- read.csv(betas_csv, stringsAsFactors = FALSE)
betas_short   <- betas[betas$RISK %in% sub(".*g__", "g__", colnames(genera_clr)), ]
betas_short$estimate <- as.numeric(betas_short$B)
genera_selection <- genera_clr[, match(betas_short$RISK,
                                       sub(".*g__", "g__", colnames(genera_clr)))]
cat("Number of genera in data:", nrow(betas_short), "\n")

risk_prediction <- genera_selection %*% betas_short$estimate
metadata2       <- subset(metadata, sampleid %in% rownames(risk_prediction))
cat("Number of participants with risk score:", nrow(metadata2), "\n")

# ============================================================
#  COX PH  (age timescale, multiple exposure scalings)
# ============================================================

if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
if (!requireNamespace("broom",    quietly = TRUE)) install.packages("broom")
if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
library(survival); library(broom); library(dplyr)

metadata2$age_entry  <- metadata2$age
metadata2$age_exit   <- metadata2$age + metadata2$studytime
metadata2$genus_FI   <- as.numeric(risk_prediction[, 1])
metadata2            <- subset(metadata2, age_exit > age_entry & !is.na(mortality))
metadata2$genus_FI_z <- as.numeric(scale(metadata2$genus_FI))

metadata2[c("sex","ppump","statin","metfor")] <-
  lapply(metadata2[c("sex","ppump","statin","metfor")], as.factor)

# ---- Null model archive (global scope via <<-) ----
null_models_archive <- list()

# ---- Helper: build per-model data subset + derived exposure ----
make_subset_with_exposure <- function(base_df, covars, exposure_kind) {
  req  <- c("age_entry","age_exit","mortality","age","genus_FI", covars)
  keep <- intersect(req, names(base_df))
  sub  <- base_df[, keep, drop = FALSE]
  sub  <- sub[complete.cases(sub), , drop = FALSE]

  fcts <- intersect(c("sex","ppump","statin","metfor"), colnames(sub))
  if (length(fcts)) sub[fcts] <- lapply(sub[fcts], as.factor)

  if      (exposure_kind == "raw")  { varname <- "genus_FI" }
  else if (exposure_kind == "z")    { varname <- "genus_FI_z"
                                      sub[[varname]] <- as.numeric(scale(sub$genus_FI)) }
  else if (exposure_kind == "AA")   { varname <- "AAgenus_FI"
                                      sub[[varname]] <- resid(lm(genus_FI ~ age, data = sub)) }
  else if (exposure_kind == "AA_z") { varname <- "AAgenus_FI_z"
                                      sub[[varname]] <- as.numeric(scale(resid(lm(genus_FI ~ age, data = sub)))) }
  else stop("Unknown exposure_kind")

  list(data = sub, varname = varname)
}

# ---- Fit one model, return tidy results + diagnostics ----
fit_one <- function(base_df, covars, exposure_kind, model_label) {
  made <- make_subset_with_exposure(base_df, covars, exposure_kind)
  df   <- made$data
  var  <- made$varname
  if (nrow(df) == 0) return(NULL)

  surv_obj <- with(df, Surv(time = age_entry, time2 = age_exit, event = mortality))

  # 1. Null model (covariates only, no biomarker)
  fmla_null <- as.formula(paste("surv_obj ~", paste(covars, collapse = " + ")))
  fit_null  <- coxph(fmla_null, data = df)
  null_models_archive[[model_label]] <<- fit_null   # archive globally

  # 2. Full model (biomarker + covariates)
  fmla_full <- as.formula(paste("surv_obj ~", paste(c(var, covars), collapse = " + ")))
  fit_full  <- coxph(fmla_full, data = df)

  # 3. Diagnostics
  lrt_p     <- anova(fit_null, fit_full)$`Pr(>|Chi|)`[2]
  zph_obj   <- cox.zph(fit_full)
  shapiro_p <- shapiro.test(zph_obj$y[, var])$p.value

  # 4. Tidy output  — N and Ncases are now both included
  out <- cbind(
    broom::tidy(fit_full, exponentiate = TRUE, conf.int = TRUE),
    model            = model_label,
    exposure         = var,
    N                = fit_full$n,          # total participants in model
    Ncases           = fit_full$nevent,     # number of events (deaths)
    concordance_null = unname(summary(fit_null)$concordance[1]),
    concordance_full = unname(summary(fit_full)$concordance[1]),
    c_index_diff     = unname(summary(fit_full)$concordance[1] -
                               summary(fit_null)$concordance[1]),
    lrt_p            = lrt_p,
    resid_norm_p     = shapiro_p
  )
  attr(out, "zph") <- zph_obj
  out
}

# ---- Run all exposure scalings across the three covariate models ----
fit_models_by_exposure <- function(base_df, exposure_kind) {
  list(
    fit_one(base_df, c("sex"),
            exposure_kind, paste0(exposure_kind, "_M1: sex")),
    fit_one(base_df, c("sex","ppump","statin","metfor"),
            exposure_kind, paste0(exposure_kind, "_M2: sex+ppump+statin+metfor")),
    fit_one(base_df, c("sex","ppump","statin","metfor","bmi"),
            exposure_kind, paste0(exposure_kind, "_M3: sex+ppump+statin+metfor+bmi"))
  )
}

print("Running Cox PH models...")
res_unit   <- fit_models_by_exposure(metadata2, "raw")
res_sd     <- fit_models_by_exposure(metadata2, "z")
res_unitAA <- fit_models_by_exposure(metadata2, "AA")
res_sdAA   <- fit_models_by_exposure(metadata2, "AA_z")
res_all    <- c(res_unit, res_sd, res_unitAA, res_sdAA)

cox_all  <- do.call(rbind, res_all)
cox_exp  <- subset(cox_all,
                   term %in% c("genus_FI","genus_FI_z","AAgenus_FI","AAgenus_FI_z"))
names(cox_exp)[names(cox_exp) %in% c("estimate","conf.low","conf.high")] <-
  c("HR","LL","UL")

# ============================================================
#  OUTPUTS
# ============================================================

# Save null model archive
saveRDS(null_models_archive,
        paste0(cohort, "_", type, "_NullModels_Archive.rds"))

# Schoenfeld residual plots
print("Saving residuals...")
pdf(paste0(cohort, "_", type, "_Residuals.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
for (i in seq_along(res_all)) {
  res_obj <- res_all[[i]]
  if (is.null(res_obj)) next
  z           <- attr(res_obj, "zph")
  var_to_plot <- res_obj$exposure[1]
  model_name  <- res_obj$model[1]
  if (!is.null(z) && var_to_plot %in% rownames(z$table)) {
    tryCatch({
      var_idx <- which(rownames(z$table) == var_to_plot)
      plot(z[var_idx], main = paste("PH Check:", model_name))
      abline(h = 0, col = "red", lty = 2)
      hist(z$y[, var_idx],
           main = paste("Residuals:", var_to_plot),
           xlab = "Schoenfeld Residuals", col = "grey")
    }, error = function(e)
      message(paste("Could not plot residuals for", model_name, ":", e$message)))
  }
}
dev.off()

write.csv(cox_exp,
          paste0(cohort, "_", type, "_CoxPH_Comprehensive.csv"),
          row.names = FALSE)

cat(sprintf(
  "Finished.\n  Results        : %s_%s_CoxPH_Comprehensive.csv\n  Residuals PDF  : %s_%s_Residuals.pdf\n  Null models    : %s_%s_NullModels_Archive.rds\n",
  cohort, type, cohort, type, cohort, type
))
