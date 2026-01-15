# Internal helper to build QIIME-style taxonomy strings
# Not exported

build_taxonomy_table <- function(tx, feature_ids) {

  pick_rank <- function(df, candidates, n) {
    col <- candidates[candidates %in% colnames(df)][1]
    if (!is.na(col)) {
      v <- as.character(df[[col]])
      v[is.na(v)] <- ""
      sub("^[a-z]__", "", v)
    } else {
      rep("", n)
    }
  }

  pref <- function(x, p) paste0(p, x)
  n <- nrow(tx)

  tax <- paste(
    pref(pick_rank(tx, c("Kingdom", "Domain"), n), "k__"),
    pref(pick_rank(tx, c("Phylum"), n), "p__"),
    pref(pick_rank(tx, c("Class"), n), "c__"),
    pref(pick_rank(tx, c("Order"), n), "o__"),
    pref(pick_rank(tx, c("Family"), n), "f__"),
    pref(pick_rank(tx, c("Genus"), n), "g__"),
    pref(pick_rank(tx, c("Species"), n), "s__"),
    sep = "; "
  )

  data.frame(
    FeatureID = feature_ids,
    Taxon = tax,
    stringsAsFactors = FALSE
  )
}

