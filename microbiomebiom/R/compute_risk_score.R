#' Compute microbiome risk score
#'
#' @description
#' Computes a genus-based microbiome risk score.
#' Users are expected to apply all cohort-specific inclusion and exclusion
#' criteria prior to analysis. Metadata is used only to align sample identifiers.
#'
#' @param input Object or path. One of:
#'   \itemize{
#'     \item QIIME2 feature table (.qza)
#'     \item phyloseq object or .rds
#'     \item TreeSummarizedExperiment object or .rds
#'   }
#' @param taxonomy_qza Path to taxonomy.qza (required if input_type = "qza")
#' @param metadata Optional data.frame with sample identifiers
#' @param sampleid_col Column name in metadata containing sample IDs
#' @param input_type One of c("qza", "phyloseq", "tse")
#' @param betas Data.frame with columns RISK and B (default: betas_default)
#'
#' @return An object of class \code{microbiomebiom_risk}
#' @export
compute_risk_score <- function(
    input,
    taxonomy_qza = NULL,
    metadata = NULL,
    sampleid_col = "sampleid",
    input_type = c("qza", "phyloseq", "tse"),
    betas = betas_default
) {

  input_type <- match.arg(input_type)

  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required")
  }

  ## ---- LOAD COUNTS + TAXONOMY ----
  if (input_type == "qza") {
    if (is.null(taxonomy_qza)) stop("taxonomy_qza must be provided for qza input")
    if (!requireNamespace("rbiom", quietly = TRUE)) stop("Package 'rbiom' is required")
    tmp <- tempfile()
    dir.create(tmp)
    system(paste("unzip -d", file.path(tmp, "counts"), input))
    system(paste("unzip -d", file.path(tmp, "taxonomy"), taxonomy_qza))
    taxonomy_table <- read.table(
      system(paste("find", file.path(tmp, "taxonomy"), "-name taxonomy.tsv"), intern = TRUE),
      header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
    biom <- rbiom::as_rbiom(
      system(paste("find", file.path(tmp, "counts"), "-name feature-table.biom"), intern = TRUE)
    )
    counts <- as.matrix(biom$counts)

  } else if (input_type == "phyloseq") {
    if (!requireNamespace("phyloseq", quietly = TRUE)) stop("phyloseq required")
    ps <- if (is.character(input)) readRDS(input) else input
    counts <- as.matrix(phyloseq::otu_table(ps))
    if (!phyloseq::taxa_are_rows(ps)) counts <- t(counts)
    tx <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
    taxonomy_table <- build_taxonomy_table(tx, rownames(tx))

  } else if (input_type == "tse") {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment required")
    tse <- if (is.character(input)) readRDS(input) else input
    counts <- SummarizedExperiment::assay(tse, 1)
    tx <- as.data.frame(SummarizedExperiment::rowData(tse), stringsAsFactors = FALSE)
    taxonomy_table <- build_taxonomy_table(tx, rownames(tse))
  }

  ## ---- SAMPLE SELECTION ----
  if (!is.null(metadata)) {
    names(metadata) <- tolower(names(metadata))
    sampleid_col <- tolower(sampleid_col)
    if (!sampleid_col %in% names(metadata)) stop("metadata does not contain column: ", sampleid_col)
    sample_ids <- as.character(metadata[[sampleid_col]])
  } else {
    sample_ids <- colnames(counts)
  }

  ## ---- COLLAPSE TO GENUS ----
  taxonomy_table$genus <- sub("; s__.*", "", taxonomy_table$Taxon)
  counts <- counts[taxonomy_table$FeatureID, , drop = FALSE]
  counts_genera <- aggregate(counts ~ factor(taxonomy_table$genus), FUN = sum)
  rownames(counts_genera) <- counts_genera[,1]
  counts_genera <- t(counts_genera[, -1])
  counts_genera <- counts_genera[intersect(sample_ids, rownames(counts_genera)), , drop = FALSE]

  ## ---- FILTER + CLR ----
  abundance <- sweep(counts_genera, 1, rowSums(counts_genera), "/")
  keep <- colSums(abundance > 0) > 0.1 * nrow(abundance) &
    apply(abundance, 2, function(x) mean(x[x > 0]) >= 1e-4)
  genera_clr <- vegan::decostand(counts_genera[, keep, drop = FALSE], method = "clr", pseudocount = 1)

  ## ---- RISK SCORE ----
  if (!all(c("RISK","B") %in% colnames(betas))) stop("betas must contain columns RISK and B")

  betas_short <- betas[betas$RISK %in% sub(".*g__", "g__", colnames(genera_clr)), ]
  betas_short$estimate <- as.numeric(betas_short$B)
  genera_sel <- genera_clr[, match(betas_short$RISK, sub(".*g__", "g__", colnames(genera_clr)))]

  risk <- as.numeric(genera_sel %*% betas_short$estimate)
  names(risk) <- rownames(genera_sel)

  ## ---- RETURN DIFFERENT OUTPUTS ----
  if (identical(betas, betas_default)) {
    out <- data.frame(
      sampleid = rownames(genera_sel),
      genusFI = risk,
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      sampleid = rownames(genera_sel),
      risk_score = risk,
      stringsAsFactors = FALSE
    )
  }

  return(out)
}
