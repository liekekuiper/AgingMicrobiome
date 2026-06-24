# ════════════════════════════════════════════════════════════════════════════
# DCS-only. Produces the two DCS PNGs that get embedded in
# the combined patchwork figure on the other server.
#
# The panels are built with the SAME styling as the live patchwork panels and
# saved at the size the patchwork script computed (dcs_panel_target.rds).
# Neither PNG carries a title: in the patchwork the cohort name is a spanning
# strip and the correlation panel has no title. The correlation PNG DOES bake
# in the R/R^2/P annotation, since the patchwork embeds it as a flat image.
#
# Workflow:
#   1. Run the patchwork script once -> it writes dcs_panel_target.rds.
#   2. Copy dcs_panel_target.rds into the wd below (this script reads it).
#   3. Run this script on RIVM server -> writes the two PNGs to the wd.
#   4. Download the PNGs, drop them into the DCS/ folder next to the patchwork
#      script, and re-run the patchwork script.
# ════════════════════════════════════════════════════════════════════════════

#### Pu et al. frailty direction per target species ##########################
frailty_direction = data.frame(
  target = c(
    "Bacteroides eggerthii", "Clostridium bolteae", "Clostridium clostridioforme",
    "Clostridium innocuum", "Erysipelatoclostridium ramosum", "Flavonifractor plautii",
    "Hungatella hathewayi", "Megamonas funiformis", "Megamonas hypermegale",
    "Parabacteroides merdae", "Prevotella sp. 885", "Proteobacteria bacterium CAG:139",
    "Ruminococcus gnavus",
    "Adlercreutzia equolifaciens", "Faecalibacterium prausnitzii",
    "Oscillibacter sp. 57 20", "Roseburia sp. CAG:471", "Ruminococcus callidus"
  ),
  direction = c(rep("positive", 13), rep("negative", 5)),
  stringsAsFactors = FALSE
)

#### Synonym table: MetaPhlAn3 -> GTDB/GG2 ##################################─
synonyms = data.frame(
  metaphlan3 = c(
    "Clostridium bolteae", "Clostridium clostridioforme", "Clostridium innocuum",
    "Ruminococcus gnavus", "Ruminococcus callidus", "Hungatella hathewayi",
    "Faecalibacterium prausnitzii"
  ),
  gg2_canonical = c(
    "Enterocloster bolteae", "Enterocloster clostridioformis", "Clostridium_AQ innocuum",
    "Ruminococcus_B gnavus", "Ruminococcus_C_58660 callidus", "Hungatella_A_127239 hathewayi_A",
    "Faecalibacterium prausnitzii_C_71351"
  ),
  stringsAsFactors = FALSE
)

#### Target species (MetaPhlAn3 nomenclature) ################################
target_species = c(
  "Bacteroides eggerthii", "Clostridium bolteae", "Clostridium clostridioforme",
   "Clostridium innocuum", "Erysipelatoclostridium ramosum",
  "Flavonifractor plautii", "Hungatella hathewayi", "Megamonas funiformis",
  "Megamonas hypermegale", "Parabacteroides merdae", "Prevotella sp. 885",
  "Proteobacteria bacterium CAG:139", "Ruminococcus gnavus",
  "Adlercreutzia equolifaciens", "Faecalibacterium prausnitzii",
  "Oscillibacter sp. 57 20", "Roseburia sp. CAG:471", "Ruminococcus callidus"
)

#### Shared palette + panel builders (identical styling to the patchwork) ####─
pal = c(No = "#7FB9A8", Mild = "#E8A0BF",
        Moderate = "#C0504D", Severe = "#6B4E8E")

fmt_p = function(p) {
  if (is.na(p)) return("P: NA")
  e = format(p, scientific = TRUE, digits = 5)
  parts = strsplit(e, "e")[[1]]
  sprintf("P: %s %%*%% 10^%d", parts[1], as.integer(parts[2]))
}

# Violin panel, NO title (cohort name is a strip in the patchwork).
make_violin_dcs = function(vdat, panel_ar = 0.85) {
  vdat$group = factor(vdat$group, levels = c("No", "Mild", "Moderate", "Severe"))
  group_n = table(vdat$group)
  present = names(group_n)[group_n > 0]
  comps = list(c("No", "Mild"), c("No", "Moderate"), c("No", "Severe"))
  comps = comps[sapply(comps, function(p) all(p %in% present))]

  labs_n = setNames(sprintf("%s\nn=%d", names(group_n), as.integer(group_n)),
                    names(group_n))

  g = ggplot(vdat, aes(x = group, y = fi_microbiota, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.85, color = NA, scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 alpha = 0.9, color = "grey25") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = labs_n)

  if (length(comps) > 0) {
    g = g + ggsignif::geom_signif(
      comparisons = comps, test = "wilcox.test",
      test.args = list(alternative = "two.sided"),
      map_signif_level = FALSE, step_increase = 0.10,
      textsize = 2.6, tip_length = 0.01)
  }

  g + labs(x = NULL, y = "Pu et al. FI-microbiota") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(color = "black"),
          aspect.ratio = panel_ar)
}

# Correlation panel, R/R^2/P annotated, NO title.
make_corr_dcs = function(cdat, R, R2, P, panel_ar = 0.80) {
  rng_x = range(cdat$fi_microbiota, na.rm = TRUE)
  rng_y = range(cdat$real_FI, na.rm = TRUE)
  ytop  = rng_y[2]; ystep = 0.10 * diff(rng_y)

  ggplot(cdat, aes(x = fi_microbiota, y = real_FI)) +
    geom_point(color = "coral", alpha = 0.6, size = 0.8) +
    geom_smooth(method = "lm", color = "#3B5BA5", fill = "grey70") +
    annotate("text", x = rng_x[1], y = ytop,
             label = sprintf("R: %s", signif(R, 5)),
             hjust = 0, vjust = 1, size = 3) +
    annotate("text", x = rng_x[1], y = ytop - ystep,
             label = sprintf("R^2: %s", signif(R2, 5)),
             hjust = 0, vjust = 1, size = 3, parse = TRUE) +
    annotate("text", x = rng_x[1], y = ytop - 2 * ystep,
             label = fmt_p(P),
             hjust = 0, vjust = 1, size = 3, parse = TRUE) +
    labs(x = "Pu et al. FI-microbiota", y = "FI") +
    theme_classic(base_size = 11) +
    theme(aspect.ratio = panel_ar)
}

#### Read the panel-size target written by the patchwork script ##############
DCS_TARGET_FILE = "dcs_panel_target.rds"

read_dcs_target = function(path) {
  default = list(
    violin = list(width_in = 3.6, height_in = 3.7, panel_ar = 0.85),
    corr   = list(width_in = 4.3, height_in = 3.7, panel_ar = 0.80),
    dpi    = 600
  )
  if (!file.exists(path)) {
    message("WARNING: ", path, " not found. Using fallback panel sizes. ",
            "Copy dcs_panel_target.rds from the patchwork server for an exact match.")
    return(default)
  }
  t = readRDS(path)
  message(sprintf("Loaded DCS target:\n  violin %.2f x %.2f in (panel_ar %.2f)\n  corr   %.2f x %.2f in (panel_ar %.2f)\n  dpi %d",
                  t$violin$width_in, t$violin$height_in, t$violin$panel_ar,
                  t$corr$width_in,   t$corr$height_in,   t$corr$panel_ar,
                  t$dpi))
  t
}

#### Main analysis function ################################################─
run_analysis_fi_microbiota <- function(wd, metafilepath, cohort, type,
                                       target_species, synonyms, frailty_direction,
                                       dcs_target) {
  setwd(wd)
  library(ggplot2)

  if (!require("rbiom")) install.packages("rbiom")
  if (!require("vegan")) install.packages("vegan")
  if (!require("ggsignif")) install.packages("ggsignif")

  #### Load data ##############################################################
  taxonomy_table = read.table(system("find ./QZA_temp/taxonomy -name taxonomy.tsv", intern = T),
                              header = T, sep = "\t")
  biom_table = rbiom::as_rbiom(system("find ./QZA_temp/counts -name feature-table.biom", intern = T))
  counts = as.matrix(biom_table$counts)

  metadata = read.table(metafilepath, header = T, sep = "\t")
  names(metadata) = tolower(names(metadata))
  metadata = metadata[!is.na(metadata$fi) & !is.na(metadata$sampleid), ]

  #### Species labels from GG2 taxonomy strings ##############################─
  extract_species_gg2 <- function(taxstring) {
    m <- regmatches(taxstring, regexpr("s__[^;]+", taxstring))
    if (length(m) == 0 || m == "s__") return(NA_character_)
    sub("s__", "", m)
  }
  taxonomy_table$species <- sapply(taxonomy_table[,2], extract_species_gg2)

  has_species <- !is.na(taxonomy_table$species) & taxonomy_table$species != ""
  taxonomy_table_species <- taxonomy_table[has_species, ]
  counts_species_all <- counts[taxonomy_table_species[,1], ]

  counts_agg_all = aggregate(counts_species_all ~ factor(taxonomy_table_species$species), FUN = sum)
  rownames(counts_agg_all) = counts_agg_all[,1]
  counts_agg_all = counts_agg_all[,-1]
  counts_agg_all = t(counts_agg_all)

  sample_ids <- as.character(metadata$sampleid)
  counts_agg_all <- counts_agg_all[sample_ids[sample_ids %in% rownames(counts_agg_all)], , drop = FALSE]

  #### Prevalence + abundance filter ##########################################
  abundance_species = sweep(counts_agg_all, 1, rowSums(counts_agg_all), '/')
  filter = colSums(abundance_species > 0) > 0.1 * nrow(abundance_species) &
    apply(abundance_species, 2, \(x) mean(x[x > 0]) >= 0.0001)
  counts_species_filtered = counts_agg_all[, filter]

  gg2_filtered = colnames(counts_species_filtered)
  gg2_all      = colnames(counts_agg_all)

  #### Match target species to GG2 ############################################
  match_results = lapply(target_species, function(tgt) {
    syn_idx  = which(synonyms$metaphlan3 == tgt)
    gg2_name = if (length(syn_idx) > 0) synonyms$gg2_canonical[syn_idx] else tgt
    in_filtered   = gg2_name %in% gg2_filtered
    in_unfiltered = gg2_name %in% gg2_all
    is_synonym_of = if (length(syn_idx) > 0 && tgt != gg2_name) {
      other = target_species[target_species != tgt &
                               (target_species == gg2_name |
                                  sapply(target_species, function(x) {
                                    si = which(synonyms$metaphlan3 == x)
                                    length(si) > 0 && synonyms$gg2_canonical[si] == gg2_name
                                  }))]
      if (length(other) > 0) paste(other, collapse = "; ") else NA
    } else NA
    list(target = tgt, gg2_match = gg2_name, synonym_used = length(syn_idx) > 0,
         synonym_of = is_synonym_of, in_filtered = in_filtered,
         in_unfiltered = in_unfiltered,
         status = ifelse(in_filtered, "in data (filtered)",
                         ifelse(in_unfiltered, "present but filtered out",
                                "absent from database")))
  })

  match_df = data.frame(
    target        = sapply(match_results, `[[`, "target"),
    gg2_match     = sapply(match_results, `[[`, "gg2_match"),
    synonym_used  = sapply(match_results, `[[`, "synonym_used"),
    synonym_of    = sapply(match_results, `[[`, "synonym_of"),
    in_filtered   = sapply(match_results, `[[`, "in_filtered"),
    in_unfiltered = sapply(match_results, `[[`, "in_unfiltered"),
    status        = sapply(match_results, `[[`, "status"),
    stringsAsFactors = FALSE
  )

  message(sprintf("\n## Species match summary for %s %s ##", cohort, type))
  message(sprintf("Targets: %d | In filtered data: %d | Filtered out: %d | Absent from DB: %d",
                  nrow(match_df), sum(match_df$in_filtered),
                  sum(!match_df$in_filtered & match_df$in_unfiltered),
                  sum(!match_df$in_unfiltered)))
  print(match_df[, c("target", "gg2_match", "synonym_of", "status")])
  write.csv(match_df, paste0(cohort, '_', type, "_pufi_species_matchdf.csv"), row.names = FALSE)

  #### FI-microbiota score ####################################################─
  match_table = merge(match_df, frailty_direction, by = "target", all.x = TRUE)
  available   = match_table[match_table$in_filtered, ]
  available   = available[!duplicated(available$gg2_match), ]

  message(sprintf("\nSpecies included in FI: %d (positive: %d, negative: %d)",
                  nrow(available),
                  sum(available$direction == "positive", na.rm = TRUE),
                  sum(available$direction == "negative", na.rm = TRUE)))

  abundance = sweep(counts_species_filtered, 1, rowSums(counts_species_filtered), '/')
  deficit_matrix = matrix(NA, nrow = nrow(abundance), ncol = nrow(available),
                          dimnames = list(rownames(abundance), available$gg2_match))
  for (i in seq_len(nrow(available))) {
    sp  = available$gg2_match[i]; dir = available$direction[i]
    if (!sp %in% colnames(abundance)) next
    x = abundance[, sp]; p0.5 = quantile(x, 0.005); p99.5 = quantile(x, 0.995)
    score = if (dir == "positive") (x - p0.5) / (p99.5 - p0.5) else (p99.5 - x) / (p99.5 - p0.5)
    deficit_matrix[, sp] = pmax(0, pmin(1, score))
  }
  fi_score = rowMeans(deficit_matrix, na.rm = TRUE)

  #### Correlation with measured FI ##########################################─
  metadata2 = subset(metadata, sampleid %in% names(fi_score))
  fi_aligned = fi_score[match(metadata2$sampleid, names(fi_score))]
  result = cor.test(fi_aligned, metadata2$fi)
  print(paste0(cohort, " | ", type, "  R:", signif(result$estimate, 5),
               "  R2:", signif(result$estimate^2, 5), "  P:", signif(result$p.value, 5)))

  cor_summary = data.frame(cohort = cohort, type = type,
                           R = signif(result$estimate, 5), R2 = signif(result$estimate^2, 5),
                           P = signif(result$p.value, 5), n = nrow(metadata2),
                           n_species = nrow(available), row.names = NULL)

  #### Build violin_df ########################################################─
  violin_df = data.frame(sampleid = metadata2$sampleid, fi_microbiota = fi_aligned,
                         fi = metadata2$fi, stringsAsFactors = FALSE)
  violin_df = violin_df[!is.na(violin_df$fi) & !is.na(violin_df$fi_microbiota), ]
  violin_df$group = cut(violin_df$fi, breaks = c(-Inf, 0.1, 0.2, 0.3, Inf),
                        labels = c("No", "Mild", "Moderate", "Severe"))
  violin_df$group = factor(violin_df$group, levels = c("No", "Mild", "Moderate", "Severe"))

  group_n = table(violin_df$group)
  message(sprintf("\n%s %s frailty group sizes: No=%d Mild=%d Moderate=%d Severe=%d",
                  cohort, type, group_n["No"], group_n["Mild"],
                  group_n["Moderate"], group_n["Severe"]))

  #### Build the two panels with patchwork-matched styling ####################─
  g_violin = make_violin_dcs(violin_df, panel_ar = dcs_target$violin$panel_ar)

  cdat = data.frame(fi_microbiota = violin_df$fi_microbiota, real_FI = violin_df$fi)
  g_corr = make_corr_dcs(cdat,
                         unname(result$estimate), unname(result$estimate^2),
                         result$p.value, panel_ar = dcs_target$corr$panel_ar)

  #### Save at the exact target panel size ####################################─
  ggsave(filename = paste0(cohort, "_", type, "_fi_microbiota_violin.png"),
         plot = g_violin,
         width = dcs_target$violin$width_in, height = dcs_target$violin$height_in,
         dpi = dcs_target$dpi)
  ggsave(filename = paste0(cohort, "_", type, "_pufi_microbiota.png"),
         plot = g_corr,
         width = dcs_target$corr$width_in, height = dcs_target$corr$height_in,
         dpi = dcs_target$dpi)

  # Deficit-distribution QC plot (not used by the patchwork).
  deficit_long = reshape(as.data.frame(deficit_matrix), varying = colnames(deficit_matrix),
                         v.names = "deficit", timevar = "species",
                         times = colnames(deficit_matrix), direction = "long")
  g2 = ggplot(deficit_long, aes(x = deficit, y = species)) +
    ggridges::geom_density_ridges(fill = "coral", alpha = 0.7) +
    labs(x = "Deficit score", y = NULL,
         title = paste0(cohort, " ", type, " — deficit score distribution per species")) +
    theme_minimal()
  ggsave(filename = paste0(cohort, "_", type, "_fi_microbiota_deficits.png"), plot = g2)

  list(cor_result = result, cor_summary = cor_summary,
       violin_data = violin_df, group_n = group_n,
       violin_plot = g_violin, corr_plot = g_corr)
}

#### Cohort configuration (DCS only) ########################################─
configs_fi = list(
  DCS_16S = list(
    wd = "/home/kuiperl/my_scratch_dir/",
    metafilepath = "./microbiome_input_file.txt",
    cohort = "DCS", type = "16S"
  )
)

#### Run ####################################################################─
dcs_target = read_dcs_target(file.path(configs_fi$DCS_16S$wd, DCS_TARGET_FILE))

results_fi = list()
cor_summary = data.frame()
for (nm in names(configs_fi)) {
  cfg = configs_fi[[nm]]
  results_fi[[nm]] = run_analysis_fi_microbiota(
    wd = cfg$wd, metafilepath = cfg$metafilepath, cohort = cfg$cohort, type = cfg$type,
    target_species = target_species, synonyms = synonyms,
    frailty_direction = frailty_direction, dcs_target = dcs_target)
  cor_summary = rbind(cor_summary, results_fi[[nm]]$cor_summary)
}

print(cor_summary)
write.csv(cor_summary, "DCS_pufi_fi_correlations.csv", row.names = FALSE)

message("DCS panels written to: ", configs_fi$DCS_16S$wd,
        "\nDownload and drop into DCS/ next to the patchwork script:\n",
        "  DCS_16S_fi_microbiota_violin.png\n  DCS_16S_pufi_microbiota.png")
