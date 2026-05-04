# microbiome_utils.R
# Unified utility functions for both TSE and phyloseq backends.
# Switch between backends by passing data_source = "TSE" or "phyloseq"
# to the main process() function.

# ============================================================
#  PACKAGES  (always required)
# ============================================================
for (.pkg in c("vegan", "data.table", "dplyr", "tibble", "tidyr")) {
  if (!requireNamespace(.pkg, quietly = TRUE)) install.packages(.pkg)
  library(.pkg, character.only = TRUE)
}
rm(.pkg)

# ============================================================
#  SHARED HELPERS
# ============================================================

#--------------------------
# find_complete
#  Filter metadata to complete cases for the given model/outcome/subset
#  and apply factor encoding. Identical logic for both backends.
#--------------------------
find_complete <- function(metadata, model, subset, outcome, factors) {
  if (!is.vector(outcome)) outcome <- c(outcome)

  vars        <- if (model != "") unlist(strsplit(model, "\\+")) else character(0)
  all_columns <- unique(c(vars, outcome))
  if ("mortality" %in% outcome) all_columns <- c(all_columns, "studytime", "age")

  meta_df <- metadata %>%
    filter(!if_any(all_of(all_columns), is.na)) %>%
    filter(age >= 18)

  if (subset == "all") {
    final_df <- meta_df
  } else if (subset %in% c("men", "women")) {
    final_df <- meta_df %>% filter(sex == subset)
  } else if (startsWith(subset, "age_")) {
    age_ranges <- list(
      age_1 = c(18, 40), age_2 = c(40, 50), age_3 = c(50, 60),
      age_4 = c(60, 70), age_5 = c(70, Inf)
    )
    rng      <- age_ranges[[subset]]
    final_df <- meta_df %>% filter(age >= rng[1], age < rng[2])
  } else {
    stop(paste("Invalid subset:", subset))
  }

  base_factors  <- c("ppump", "metfor", "statin", "race")
  factor_names  <- if (length(factors) == 1 && is.na(factors)) base_factors else c(base_factors, factors)
  if (!(subset %in% c("men", "women"))) factor_names <- c("sex", factor_names)

  final_df <- final_df %>% mutate(across(all_of(factor_names), as.factor))
  return(final_df)
}

#--------------------------
# calculate_min_dissimilarity
#  Minimum pairwise dissimilarity per sample (row) from a dist object.
#--------------------------
calculate_min_dissimilarity <- function(distmat) {
  mat        <- as.matrix(distmat)
  diag(mat)  <- NA
  apply(mat, 1, function(x) min(x, na.rm = TRUE))
}

# ============================================================
#  TSE-SPECIFIC HELPERS
# ============================================================

# pick_rank: return first matching column name from candidates, or NULL
.pick_rank <- function(tx, candidates) {
  for (cand in candidates) {
    if (cand %in% colnames(tx)) return(cand)
  }
  NULL
}

# pref: package a column name + prefix into a list for aggregate_by_rank
.pref <- function(col, prefix) {
  if (is.null(col)) return(NULL)
  list(col = col, prefix = prefix)
}

#--------------------------
# aggregate_by_rank
#  Collapse a TSE assay to a given taxonomic rank with abundance filtering.
#  lenient = FALSE  →  strict:  rowMeans(rel > 0.01) >= 0.1
#  lenient = TRUE   →  lenient: prevalence > 10 % & mean positive abundance >= 1e-4
#--------------------------
aggregate_by_rank <- function(tse_obj, assay_name = "counts", rank_info, lenient = FALSE) {
  if (is.null(rank_info)) stop("No suitable taxonomy column found")

  x     <- SummarizedExperiment::assays(tse_obj)[[assay_name]]
  tax   <- SummarizedExperiment::rowData(tse_obj)
  group <- as.character(tax[[rank_info$col]])
  group[is.na(group) | group == ""] <- "unclassified"

  agg       <- rowsum(x, group = group, reorder = FALSE)
  new_names <- rownames(agg)
  if (!is.null(rank_info$prefix)) new_names <- paste0(rank_info$prefix, new_names)
  rownames(agg) <- new_names

  lib <- colSums(agg); lib[lib == 0] <- 1
  rel <- sweep(agg, 2, lib, "/")

  if (lenient) {
    keep <- rowSums(rel > 0) > 0.1 * ncol(rel) &
      apply(rel, 1, function(v) {
        vp <- v[v > 0]; if (length(vp) == 0) FALSE else mean(vp) >= 1e-4
      })
  } else {
    keep <- rowMeans(rel > 0.01) >= 0.1
  }

  agg[keep, , drop = FALSE]
}

#--------------------------
# clr_transform  (TSE path)
#  Manual CLR on a features-x-samples count matrix.
#--------------------------
clr_transform <- function(counts_mat, pseudocount = 1) {
  X    <- counts_mat + pseudocount
  X    <- sweep(X, 2, colSums(X), "/")
  Xt   <- t(X)
  clr  <- t(apply(Xt, 1, function(v) { lv <- log(v); lv - mean(lv) }))
  as.data.frame(clr, check.names = FALSE)
}

#--------------------------
# alpha_diversity_tse  (TSE path)
#  Thin wrapper around vegan::diversity for shannon / simpson / invsimpson.
#--------------------------
alpha_diversity_tse <- function(counts_mat, metric = c("shannon", "simpson", "invsimpson")) {
  metric <- match.arg(metric)
  vegan::diversity(t(counts_mat), index = metric)
}

# ============================================================
#  PHYLOSEQ-SPECIFIC HELPERS
# ============================================================

#--------------------------
# safe_tax_glom_filter  (phyloseq path)
#  Agglomerate at a taxonomic rank, apply abundance filter, assign prefixed names.
#  lenient = FALSE  →  strict:  rowMeans(otu_table > 0.01) >= 0.1
#  lenient = TRUE   →  lenient: prevalence > 10 % & mean positive abundance >= 1e-4
#--------------------------
safe_tax_glom_filter <- function(physeq_obj, rank, lenient = FALSE) {
  if (!(rank %in% phyloseq::rank_names(physeq_obj)))
    stop(paste("Taxonomic rank", rank, "not found in tax_table"))

  ps_glom <- phyloseq::tax_glom(physeq_obj, taxrank = rank)
  ps_rel  <- phyloseq::transform_sample_counts(ps_glom, function(x) x / sum(x))

  mat <- as(phyloseq::otu_table(ps_rel), "matrix")
  if (!phyloseq::taxa_are_rows(ps_rel)) mat <- t(mat)

  if (lenient) {
    keep <- rowSums(mat > 0) > 0.1 * ncol(mat) &
      apply(mat, 1, function(v) {
        vp <- v[v > 0]; if (length(vp) == 0) FALSE else mean(vp) >= 1e-4
      })
  } else {
    keep <- rowMeans(mat > 0.01) >= 0.1
  }

  if (sum(keep) == 0) stop(paste("No taxa left after filtering at", rank, "level"))

  ps_filt      <- phyloseq::prune_taxa(rownames(mat)[keep], ps_rel)
  taxa_labels  <- as.character(phyloseq::tax_table(ps_filt)[, rank])
  taxa_labels[is.na(taxa_labels)] <- "unclassified"
  prefix       <- paste0(tolower(substr(rank, 1, 1)), "__")
  phyloseq::taxa_names(ps_filt) <- paste0(prefix, taxa_labels)

  return(ps_filt)
}

#--------------------------
# to_clr  (phyloseq path)
#  CLR transform via microbiome::transform.
#--------------------------
to_clr <- function(physeq_obj) {
  clr_obj <- microbiome::transform(physeq_obj, "clr")
  as.data.frame(t(phyloseq::otu_table(clr_obj)))
}

#--------------------------
# add_alpha_to_meta  (phyloseq path)
#  Add one alpha-diversity column to a metadata data frame.
#--------------------------
add_alpha_to_meta <- function(meta_df, physeq_obj, metric, colname) {
  alpha_vals <- phyloseq::estimate_richness(physeq_obj, measures = metric)
  alpha_df   <- tibble(.rowid = rownames(alpha_vals), value = alpha_vals[[metric]])
  meta_df %>%
    rownames_to_column(".rowid") %>%
    left_join(alpha_df, by = ".rowid") %>%
    rename(!!colname := value) %>%
    column_to_rownames(".rowid")
}

# ============================================================
#  UNIFIED PROCESS FUNCTION
# ============================================================

#' process
#'
#' @param data_obj    A TSE or phyloseq object (must match data_source).
#' @param metadata_path  Path to tab-delimited metadata file.
#' @param model       Covariate string, e.g. "sex+age".
#' @param sub         Subset label: "all", "men", "women", or "age_1"–"age_5".
#' @param outcome     Outcome variable name, e.g. "mortality".
#' @param factors     Character vector of cohort-specific factor columns.
#' @param label       "16s" (genus only) or other string (adds species level).
#' @param data_source "TSE" or "phyloseq".
#' @param assay_name  Assay name in TSE (ignored for phyloseq). Default "counts".
#' @param lenient     Logical. Use lenient abundance filter? Default FALSE.
#'
#' @return A data frame combining metadata with CLR features and diversity metrics.
process <- function(data_obj,
                    metadata_path,
                    model,
                    sub,
                    outcome,
                    factors,
                    label       = "16s",
                    data_source,
                    assay_name  = "counts",
                    lenient     = FALSE) {

  if (!data_source %in% c("TSE", "phyloseq"))
    stop(sprintf("Unknown data_source '%s'. Use \"TSE\" or \"phyloseq\".", data_source))

  # ---- Load & filter metadata (shared) ----
  meta <- data.table::fread(metadata_path) %>% as.data.frame()
  colnames(meta)  <- tolower(colnames(meta))
  meta$sampleid   <- as.character(meta$sampleid)
  rownames(meta)  <- meta$sampleid

  meta_df <- find_complete(meta, model, sub, outcome, factors)

  # ===========================================================
  if (data_source == "TSE") {
  # ===========================================================

    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
      BiocManager::install("SummarizedExperiment")
    library(SummarizedExperiment)

    # Subset TSE to common sample IDs
    counts      <- SummarizedExperiment::assays(data_obj)[[assay_name]]
    common_ids  <- intersect(colnames(counts), rownames(meta_df))
    tse_sub     <- data_obj[, common_ids]
    counts      <- SummarizedExperiment::assays(tse_sub)[[assay_name]]
    meta_df     <- meta_df[common_ids, , drop = FALSE]

    tx    <- SummarizedExperiment::rowData(tse_sub)
    tax_G <- .pref(.pick_rank(tx, c("Genus","genus","Rank6")),   "g__")
    tax_S <- .pref(.pick_rank(tx, c("Species","species","Rank7")), "s__")

    # Genus CLR
    genus_counts <- aggregate_by_rank(tse_sub, assay_name, tax_G, lenient)
    genus_clr    <- clr_transform(genus_counts)
    genus_clr$.rowid <- rownames(genus_clr)

    # Species CLR (WGS only)
    species_counts <- NULL
    species_clr    <- NULL
    if (tolower(label) != "16s" && !is.null(tax_S)) {
      species_counts <- aggregate_by_rank(tse_sub, assay_name, tax_S, lenient)
      species_clr    <- clr_transform(species_counts)
      species_clr$.rowid <- rownames(species_clr)
    }

    # Alpha diversity
    alpha_metrics <- c("shannon", "simpson", "invsimpson")
    for (m in alpha_metrics) {
      meta_df[[paste0(m, "_asv")]]   <- alpha_diversity_tse(counts,        m)[rownames(meta_df)]
      meta_df[[paste0(m, "_genus")]] <- alpha_diversity_tse(genus_counts,  m)[rownames(meta_df)]
      if (!is.null(species_counts))
        meta_df[[paste0(m, "_species")]] <- alpha_diversity_tse(species_counts, m)[rownames(meta_df)]
    }

    # Beta diversity
    meta_df$min_bray_asv   <- calculate_min_dissimilarity(vegan::vegdist(t(counts),        "bray"))
    meta_df$min_jacc_asv   <- calculate_min_dissimilarity(vegan::vegdist(t(counts),        "jaccard"))
    meta_df$min_bray_genus <- calculate_min_dissimilarity(vegan::vegdist(t(genus_counts),  "bray"))
    meta_df$min_jacc_genus <- calculate_min_dissimilarity(vegan::vegdist(t(genus_counts),  "jaccard"))
    if (!is.null(species_counts)) {
      meta_df$min_bray_species <- calculate_min_dissimilarity(vegan::vegdist(t(species_counts), "bray"))
      meta_df$min_jacc_species <- calculate_min_dissimilarity(vegan::vegdist(t(species_counts), "jaccard"))
    }

    # Merge CLR features into metadata
    meta_df <- meta_df %>%
      rownames_to_column(".rowid") %>%
      left_join(genus_clr, by = ".rowid") %>%
      { if (!is.null(species_clr)) left_join(., species_clr, by = ".rowid") else . } %>%
      column_to_rownames(".rowid")

  # ===========================================================
  } else {  # phyloseq
  # ===========================================================

    if (!requireNamespace("phyloseq",    quietly = TRUE)) install.packages("phyloseq")
    if (!requireNamespace("microbiome",  quietly = TRUE)) install.packages("microbiome")
    library(phyloseq); library(microbiome)

    # Subset phyloseq to common sample IDs
    common_ids        <- intersect(phyloseq::sample_names(data_obj), rownames(meta_df))
    physeq_sub        <- phyloseq::prune_samples(common_ids, data_obj)
    meta_df           <- meta_df[common_ids, , drop = FALSE]
    phyloseq::sample_data(physeq_sub) <- phyloseq::sample_data(meta_df)

    # Genus agglomeration
    genus_physeq_counts <- phyloseq::tax_glom(physeq_sub, taxrank = "Genus")
    genus_physeq_clr    <- safe_tax_glom_filter(physeq_sub, "Genus", lenient)
    genus_clr           <- to_clr(genus_physeq_clr)

    # Species agglomeration (WGS only)
    species_physeq_counts <- NULL
    species_clr           <- NULL
    if (tolower(label) != "16s") {
      species_physeq_counts <- phyloseq::tax_glom(physeq_sub, taxrank = "Species")
      species_physeq_clr    <- safe_tax_glom_filter(physeq_sub, "Species", lenient)
      species_clr           <- to_clr(species_physeq_clr)
    }

    # Alpha diversity
    # Note: phyloseq also computes Chao1 which is not available via vegan;
    # Shannon / Simpson / InvSimpson match the TSE output column names.
    alpha_metrics <- c("Shannon", "Simpson", "InvSimpson", "Chao1")
    for (m in alpha_metrics) {
      col_suffix <- tolower(m)   # e.g. "shannon", "invsimpson", "chao1"
      meta_df <- add_alpha_to_meta(meta_df, physeq_sub,          m, paste0(col_suffix, "_asv"))
      meta_df <- add_alpha_to_meta(meta_df, genus_physeq_counts, m, paste0(col_suffix, "_genus"))
      if (!is.null(species_physeq_counts))
        meta_df <- add_alpha_to_meta(meta_df, species_physeq_counts, m, paste0(col_suffix, "_species"))
    }

    # Beta diversity
    meta_df$min_bray_asv   <- calculate_min_dissimilarity(phyloseq::distance(physeq_sub,          "bray"))
    meta_df$min_jacc_asv   <- calculate_min_dissimilarity(phyloseq::distance(physeq_sub,          "jaccard"))
    meta_df$min_bray_genus <- calculate_min_dissimilarity(phyloseq::distance(genus_physeq_counts, "bray"))
    meta_df$min_jacc_genus <- calculate_min_dissimilarity(phyloseq::distance(genus_physeq_counts, "jaccard"))
    if (!is.null(species_physeq_counts)) {
      meta_df$min_bray_species <- calculate_min_dissimilarity(phyloseq::distance(species_physeq_counts, "bray"))
      meta_df$min_jacc_species <- calculate_min_dissimilarity(phyloseq::distance(species_physeq_counts, "jaccard"))
    }

    # Merge CLR features into metadata
    genus_clr$.rowid <- rownames(genus_clr)
    meta_df <- meta_df %>%
      rownames_to_column(".rowid") %>%
      left_join(genus_clr, by = ".rowid") %>%
      column_to_rownames(".rowid")

    if (!is.null(species_clr)) {
      species_clr$.rowid <- rownames(species_clr)
      meta_df <- meta_df %>%
        rownames_to_column(".rowid") %>%
        left_join(species_clr, by = ".rowid") %>%
        column_to_rownames(".rowid")
    }
  }

  return(meta_df)
}
