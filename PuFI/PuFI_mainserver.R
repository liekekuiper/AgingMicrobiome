#### Pu et al. Frailty direction per target species ##################################
# Based on the original study results section:
# 13 species increased in abundance with frailty progression (positive association)
# 5 species decreased in abundance with frailty progression (negative association)
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

#### Synonym table: known NCBI/MetaPhlAn3 to GTDB/GG2 renames ############─
# GG2 uses GTDB taxonomy, which reclassifies many classically named genera
# into new monophyletic genera. This table maps MetaPhlAn3 names to their
# GG2 canonical equivalents so that target species can be correctly identified.
synonyms = data.frame(
  metaphlan3 = c(
    "Clostridium bolteae",
    "Clostridium clostridioforme",
    "Clostridium innocuum",
    "Ruminococcus gnavus",
    "Ruminococcus callidus",
    "Hungatella hathewayi",
    "Faecalibacterium prausnitzii"
  ),
  gg2_canonical = c(
    "Enterocloster bolteae",
    "Enterocloster clostridioformis",
    "Clostridium_AQ innocuum",
    "Ruminococcus_B gnavus",
    "Ruminococcus_C_58660 callidus",
    "Hungatella_A_127239 hathewayi_A",
    "Faecalibacterium prausnitzii_C_71351"
  ),
  stringsAsFactors = FALSE
)

#### Target species from the original study (MetaPhlAn3 nomenclature) ######
target_species = c(
  "Bacteroides eggerthii", "Clostridium bolteae", "Clostridium clostridioforme",
  "Clostridium innocuum", "Erysipelatoclostridium ramosum",
  "Flavonifractor plautii", "Hungatella hathewayi", "Megamonas funiformis",
  "Megamonas hypermegale", "Parabacteroides merdae", "Prevotella sp. 885",
  "Proteobacteria bacterium CAG:139", "Ruminococcus gnavus",
  "Adlercreutzia equolifaciens", "Faecalibacterium prausnitzii",
  "Oscillibacter sp. 57 20", "Roseburia sp. CAG:471", "Ruminococcus callidus"
)

#### Main analysis function ################################################─
run_analysis_fi_microbiota <- function(wd, metafilepath, cohort, type,
                                       target_species, synonyms, frailty_direction) {
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
  
  #### Extract species-level labels from Greengenes2 taxonomy strings ########─
  # GG2 format: "d__...; p__...; ...; g__Genus; s__Genus species"
  # Returns NA for ASVs without a species-level annotation
  extract_species_gg2 <- function(taxstring) {
    m <- regmatches(taxstring, regexpr("s__[^;]+", taxstring))
    if (length(m) == 0 || m == "s__") return(NA_character_)
    sub("s__", "", m)
  }
  
  taxonomy_table$species <- sapply(taxonomy_table[,2], extract_species_gg2)
  
  # Remove ASVs without species-level annotation
  has_species <- !is.na(taxonomy_table$species) & taxonomy_table$species != ""
  taxonomy_table_species <- taxonomy_table[has_species, ]
  counts_species_all <- counts[taxonomy_table_species[,1], ]
  
  # Aggregate ASV counts to species level
  counts_agg_all = aggregate(counts_species_all ~ factor(taxonomy_table_species$species), FUN = sum)
  rownames(counts_agg_all) = counts_agg_all[,1]
  counts_agg_all = counts_agg_all[,-1]
  counts_agg_all = t(counts_agg_all)
  
  # Retain only samples present in metadata
  sample_ids <- as.character(metadata$sampleid)
  counts_agg_all <- counts_agg_all[sample_ids[sample_ids %in% rownames(counts_agg_all)], , drop = FALSE]
  
  #### Prevalence and abundance filter ######################################─
  # Keep species present in >10% of samples with mean relative abundance >=0.01%
  # among samples where the species is detected
  abundance_species = sweep(counts_agg_all, 1, rowSums(counts_agg_all), '/')
  
  filter = colSums(abundance_species > 0) > 0.1 * nrow(abundance_species) &
    apply(abundance_species, 2, \(x) mean(x[x > 0]) >= 0.0001)
  
  counts_species_filtered = counts_agg_all[, filter]
  
  gg2_filtered = colnames(counts_species_filtered)
  gg2_all      = colnames(counts_agg_all)
  
  #### Match target species to GG2 species ##################################─
  # Apply synonym lookup first, then check availability in filtered and
  # unfiltered data. Duplicate GG2 matches across targets are allowed and
  # flagged: e.g. Clostridium bolteae and Enterocloster bolteae both map to
  # Enterocloster bolteae in GG2 and represent the same biological species.
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
    
    list(
      target        = tgt,
      gg2_match     = gg2_name,
      synonym_used  = length(syn_idx) > 0,
      synonym_of    = is_synonym_of,
      in_filtered   = in_filtered,
      in_unfiltered = in_unfiltered,
      status        = ifelse(in_filtered, "in data (filtered)",
                             ifelse(in_unfiltered, "present but filtered out",
                                    "absent from database"))
    )
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
  write.csv(match_df, paste0(cohort, '_', type, "_pufi_species_matchdf.csv"), row.names = FALSE)
  
  message(sprintf("\n## Species match summary for %s %s ##", cohort, type))
  message(sprintf("Targets: %d | In filtered data: %d | Filtered out: %d | Absent from DB: %d",
                  nrow(match_df),
                  sum(match_df$in_filtered),
                  sum(!match_df$in_filtered & match_df$in_unfiltered),
                  sum(!match_df$in_unfiltered)))
  print(match_df[, c("target", "gg2_match", "synonym_of", "status")])
  
  #### Compute FI-microbiota score ############################################
  # Implements the deficit accumulation approach from the original study:
  # - Positive association: below 0.5th percentile -> 0, above 99.5th -> 1
  # - Negative association: scoring reversed
  # - Linear interpolation between percentiles, clipped to [0, 1]
  # - Individual FI = mean deficit across all available species
  match_table = merge(match_df, frailty_direction, by = "target", all.x = TRUE)
  available   = match_table[match_table$in_filtered, ]
  available   = available[!duplicated(available$gg2_match), ]  # one entry per GG2 species
  
  message(sprintf("\nSpecies included in FI: %d (positive: %d, negative: %d)",
                  nrow(available),
                  sum(available$direction == "positive", na.rm = TRUE),
                  sum(available$direction == "negative", na.rm = TRUE)))
  
  abundance = sweep(counts_species_filtered, 1, rowSums(counts_species_filtered), '/')
  
  deficit_matrix = matrix(NA,
                          nrow = nrow(abundance),
                          ncol = nrow(available),
                          dimnames = list(rownames(abundance), available$gg2_match))
  
  for (i in seq_len(nrow(available))) {
    sp  = available$gg2_match[i]
    dir = available$direction[i]
    if (!sp %in% colnames(abundance)) next
    
    x     = abundance[, sp]
    p0.5  = quantile(x, 0.005)
    p99.5 = quantile(x, 0.995)
    
    if (dir == "positive") {
      score = (x - p0.5) / (p99.5 - p0.5)
    } else {
      score = (p99.5 - x) / (p99.5 - p0.5)
    }
    deficit_matrix[, sp] = pmax(0, pmin(1, score))
  }
  
  fi_score = rowMeans(deficit_matrix, na.rm = TRUE)
  
  #### Correlation with measured FI ##########################################─
  metadata2 = subset(metadata, sampleid %in% names(fi_score))
  fi_aligned = fi_score[match(metadata2$sampleid, names(fi_score))]
  
  result = cor.test(fi_aligned, metadata2$fi)
  print(paste0(cohort, " | ", type,
               "  R:",  signif(result$estimate, 5),
               "  R2:", signif(result$estimate^2, 5),
               "  P:",  signif(result$p.value, 5)))
  
  cor_summary = data.frame(
    cohort    = cohort,
    type      = type,
    R         = signif(result$estimate, 5),
    R2        = signif(result$estimate^2, 5),
    P         = signif(result$p.value, 5),
    n         = nrow(metadata2),
    n_species = nrow(available),
    row.names = NULL
  )
  
  #### Plot: FI-microbiota vs clinical FI ####################################
  plot_data = data.frame(fi_microbiota = fi_aligned, real_FI = metadata2$fi)
  
  g1 = ggplot(plot_data, aes(x = fi_microbiota, y = real_FI)) +
    geom_point(color = "coral") +
    geom_smooth(method = "lm") +
    labs(x = "Pu et al. FI-microbiota",
         y = "FI",
         title = paste0(cohort, " ", type, " — PuFI vs FI")) +
    theme_minimal()
  
  ggsave(filename = paste0(cohort, "_", type, "_pufi_microbiota.png"), plot = g1)
  
  #### Plot: deficit score distribution per species ##########################─
  deficit_long = reshape(
    as.data.frame(deficit_matrix),
    varying   = colnames(deficit_matrix),
    v.names   = "deficit",
    timevar   = "species",
    times     = colnames(deficit_matrix),
    direction = "long"
  )
  
  g2 = ggplot(deficit_long, aes(x = deficit, y = species)) +
    ggridges::geom_density_ridges(fill = "coral", alpha = 0.7) +
    labs(x = "Deficit score", y = NULL,
         title = paste0(cohort, " ", type, " — deficit score distribution per species")) +
    theme_minimal()
  
  ggsave(filename = paste0(cohort, "_", type, "_fi_microbiota_deficits.png"), plot = g2)
  
  #### Plot: FI-microbiota violin across frailty groups ######################─
  # Replicates Pu et al. Extended Data Fig. 10a/c. Participants are grouped by
  # measured FI using the Rockwood cut-offs (no <=0.1, mild 0.1-0.2,
  # moderate 0.2-0.3, severe >0.3). Pairwise two-sided Wilcoxon tests compare
  # each frailty group against the no-frailty reference.
  violin_df = data.frame(
    sampleid      = metadata2$sampleid,
    fi_microbiota = fi_aligned,
    fi            = metadata2$fi,
    stringsAsFactors = FALSE
  )
  violin_df = violin_df[!is.na(violin_df$fi) & !is.na(violin_df$fi_microbiota), ]
  
  violin_df$group = cut(violin_df$fi,
                        breaks = c(-Inf, 0.1, 0.2, 0.3, Inf),
                        labels = c("No", "Mild", "Moderate", "Severe"))
  violin_df$group = factor(violin_df$group,
                           levels = c("No", "Mild", "Moderate", "Severe"))
  
  group_n = table(violin_df$group)
  message(sprintf("\n%s %s frailty group sizes: No=%d Mild=%d Moderate=%d Severe=%d",
                  cohort, type,
                  group_n["No"], group_n["Mild"],
                  group_n["Moderate"], group_n["Severe"]))
  
  # Only include comparisons for groups that are actually present and non-empty
  present_groups = names(group_n)[group_n > 0]
  comps = list(c("No", "Mild"), c("No", "Moderate"), c("No", "Severe"))
  comps = comps[sapply(comps, function(p) all(p %in% present_groups))]
  
  pal = c(No = "#7FB9A8", Mild = "#E8A0BF",
          Moderate = "#C0504D", Severe = "#6B4E8E")
  
  group_labels = setNames(
    sprintf("%s\nn=%d", names(group_n), as.integer(group_n)),
    names(group_n)
  )
  g3 = ggplot(violin_df, aes(x = group, y = fi_microbiota, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.85, color = NA, scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 alpha = 0.9, color = "grey25") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = group_labels)
  
  if (length(comps) > 0) {
    g3 = g3 + ggsignif::geom_signif(
      comparisons = comps,
      test = "wilcox.test",
      test.args = list(alternative = "two.sided"),
      map_signif_level = FALSE,
      step_increase = 0.10,
      textsize = 3, tip_length = 0.01
    )
  }
  
  g3 = g3 +
    labs(x = NULL, y = "Pu et al. FI-microbiota",
         title = paste0(cohort, " ", type)) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(color = "black"))
  
  ggsave(filename = paste0(cohort, "_", type, "_fi_microbiota_violin.png"),
         plot = g3, width = 4, height = 5, dpi = 300)
  
  return(list(
    counts_filtered   = counts_species_filtered,
    counts_unfiltered = counts_agg_all,
    match_table       = match_df,
    deficit_matrix    = deficit_matrix,
    fi_score          = fi_score,
    species_used      = available,
    cor_result        = result,
    cor_summary       = cor_summary,
    all_species       = gg2_all,
    violin_data       = violin_df,
    group_n           = group_n,
    violin_plot       = g3
  ))
}

#### Cohort configurations ##################################################
configs_fi = list(
  
  MrOS_16S = list(
    wd = "~/mirna_files/Microbiome/MrOS/16S/",
    metafilepath = "../MrOS_meta.txt",
    cohort = "MrOS",
    type = "16S"
  ),
  
  MrOS_MG = list(
    wd = "~/mirna_files/Microbiome/MrOS/MrOS_Shotgun_Sequencing/",
    metafilepath = "../MrOS_meta.txt",
    cohort = "MrOS",
    type = "MG"
  ),
  FHS_16S = list(
    wd = "~/mirna_files/Microbiome/FHS/16S/",
    metafilepath = "../fhs_metadata.txt",
    cohort = "FHS",
    type = "16S"
  ),
  
  FHS_MG = list(
    wd = "~/mirna_files/Microbiome/FHS/metagenomics/",
    metafilepath = "../fhs_metadata.txt",
    cohort = "FHS",
    type = "MG"
  ),
  
  RS_16S = list(
    wd = "~/mirna_files/Microbiome/Analyses/",
    metafilepath = "metadata_agingmicrobiome.txt",
    cohort = "RS",
    type = "16S"
  )
)

#### Run across cohorts ####################################################─
results_fi = list()
cor_summary = data.frame()

for (nm in names(configs_fi)) {
  cfg = configs_fi[[nm]]
  results_fi[[nm]] = run_analysis_fi_microbiota(
    wd                = cfg$wd,
    metafilepath      = cfg$metafilepath,
    cohort            = cfg$cohort,
    type              = cfg$type,
    target_species    = target_species,
    synonyms          = synonyms,
    frailty_direction = frailty_direction
  )
  cor_summary = rbind(cor_summary, results_fi[[nm]]$cor_summary)
}

setwd('~/Microbiome/')
print(cor_summary)
write.csv(cor_summary, "pufi_fi_correlations.csv", row.names = FALSE)

DCS_violin = 'DCS/DCS_16S_fi_microbiota_violin.png'
DCS_coplot = 'DCS/DCS_16S_pufi_microbiota.png'
DCS_defict = 'DCS/DCS_16S_fi_microbiota_deficits.png'
DCS_corres = 'DCS/DCS_pufi_fi_correlations.csv'
# ════════════════════════════════════════════════════════════════════════════
# Combined figure: 6 cohorts, 2 per row, 3 rows.
# Each cohort = [violin | correlation]. Each row = [V C | V C] (4 panels),
# 12 panels total, tagged a–l. Cohort name centred above each [V C] block.
#
# DCS exists only as saved PNGs (other server). They are NOT re-plotted here;
# they are embedded as images. To make them match the live panels, the DCS
# script saves them at exactly the panel size this script computes (written to
# dcs_panel_target.rds). The sizing is fully DETERMINISTIC — no gtable measuring
# (that was unreliable). If DCS still looks a touch short/tall, nudge
# dcs_height_calib.
#
# Run AFTER the main analysis loop, so `results_fi` is in memory.
# ════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(patchwork)
if (!require("magick")) install.packages("magick")

# ── Display names (control alphabetical ordering and panel titles) ───────────
display_names = c(
  MrOS_16S = "MrOS 16S",
  MrOS_MG  = "MrOS MG",
  FHS_16S  = "Framingham Heart Study 16S",
  FHS_MG   = "Framingham Heart Study MG",
  RS_16S   = "Rotterdam Study"
)

pal = c(No = "#7FB9A8", Mild = "#E8A0BF",
        Moderate = "#C0504D", Severe = "#6B4E8E")

fmt_p = function(p) {
  if (is.na(p)) return("P: NA")
  e = format(p, scientific = TRUE, digits = 5)
  parts = strsplit(e, "e")[[1]]
  sprintf("P: %s %%*%% 10^%d", parts[1], as.integer(parts[2]))
}

# ════════════════════════════════════════════════════════════════════════════
# TUNING KNOBS
# ════════════════════════════════════════════════════════════════════════════
ar_violin = 0.95    # height/width of a violin panel. LOWER = wider.
ar_corr   = 0.80    # height/width of a correlation panel.

w_violin  = 1.0     # relative widths of the two panels within a cohort block
w_corr    = 1.2

fig_width_in  = 15  # whole-figure device size (width split across 2 cohorts)
row_height_in = 4.2
out_dpi       = 600

title_strip_h = 0.12   # block height fraction given to the cohort name strip
tag_x = 0.02           # tag (a..l) position within each panel
tag_y = 0.92

dcs_target_file = "dcs_panel_target.rds"

# Crop fractions off the DCS PNG edges before embedding. The PNG carries more
# white margin TOP/BOTTOM than left/right, so trim_frac_y is usually larger.
dcs_trim_frac_x = 0.02
dcs_trim_frac_y = 0.1

# After trimming, stretch the image to fill the FULL cell (width+height = 1 npc)
# so the DCS plot grows to match its neighbours instead of leaving a gap.
dcs_fill_cell = TRUE

# Only used when dcs_fill_cell = FALSE: 0 = flush left, 0.5 = centred.
dcs_violin_halign = 0

# ── DCS reference PNGs ──────────────────────────────────────────────────────
dcs_violin = "DCS/DCS_16S_fi_microbiota_violin.png"
dcs_coplot = "DCS/DCS_16S_pufi_microbiota.png"

# ── Re-plot a violin panel (no title; cohort name is a strip above) ──────────
make_violin = function(vdat) {
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
          aspect.ratio = ar_violin)
}

# ── Re-plot a correlation panel ──────────────────────────────────────────────
make_corr = function(cdat, R, R2, P) {
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
    theme(aspect.ratio = ar_corr)
}

# ── Embed a saved DCS PNG as a panel ─────────────────────────────────────────
# trim_frac_x / trim_frac_y crop a uniform fraction off the left+right and
# top+bottom edges respectively. The DCS PNG has more white margin top/bottom,
# so trim_frac_y is usually larger than trim_frac_x.
#
# fill_cell = TRUE stretches the (trimmed) image to the FULL cell (width and
# height = 1 npc), so after trimming the plot grows to fit the cell instead of
# keeping its own aspect ratio and leaving a gap. Because trimming removes the
# white margin (not the plot), the stretch is small and distortion minimal.
#
# halign only matters when fill_cell = FALSE: 0 = flush left, 0.5 = centred.
png_panel = function(path, trim_frac_x = 0, trim_frac_y = 0,
                     fill_cell = TRUE, halign = 0.5) {
  if (!file.exists(path)) {
    return(wrap_elements(grid::textGrob(paste0("missing:\n", path))))
  }
  img = magick::image_read(path)
  if (trim_frac_x > 0 || trim_frac_y > 0) {
    info = magick::image_info(img)
    dx = round(trim_frac_x * info$width)
    dy = round(trim_frac_y * info$height)
    geom = sprintf("%dx%d+%d+%d", info$width - 2 * dx, info$height - 2 * dy, dx, dy)
    img = magick::image_crop(img, geom)
  }
  if (fill_cell) {
    rg = grid::rasterGrob(img, interpolate = TRUE,
                          width  = grid::unit(1, "npc"),
                          height = grid::unit(1, "npc"))
  } else {
    rg = grid::rasterGrob(img, interpolate = TRUE,
                          x = grid::unit(halign, "npc"), hjust = halign)
  }
  ggplot() +
    annotation_custom(rg, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void() +
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "pt"))
}

# ── Build records for live cohorts ───────────────────────────────────────────
records = list()
for (nm in names(results_fi)) {
  r = results_fi[[nm]]
  disp = if (!is.na(display_names[nm])) display_names[nm] else nm
  cdat = data.frame(fi_microbiota = r$violin_data$fi_microbiota,
                    real_FI       = r$violin_data$fi)
  records[[disp]] = list(
    name   = disp,
    violin = make_violin(r$violin_data),
    corr   = make_corr(cdat,
                       unname(r$cor_result$estimate),
                       unname(r$cor_result$estimate^2),
                       r$cor_result$p.value)
  )
}

# ── Add DCS from saved PNGs ──────────────────────────────────────────────────
if (file.exists(dcs_violin) && file.exists(dcs_coplot)) {
  records[["Doetinchem Cohort Study"]] = list(
    name   = "Doetinchem Cohort Study",
    violin = png_panel(dcs_violin, trim_frac_x = dcs_trim_frac_x,
                       trim_frac_y = dcs_trim_frac_y,
                       fill_cell = dcs_fill_cell, halign = dcs_violin_halign),
    corr   = png_panel(dcs_coplot, trim_frac_x = dcs_trim_frac_x,
                       trim_frac_y = dcs_trim_frac_y,
                       fill_cell = dcs_fill_cell)
  )
} else {
  message("DCS PNGs not found — DCS row skipped. Expected: ",
          dcs_violin, " and ", dcs_coplot)
}

# ── Order alphabetically by display name ────────────────────────────────────
records = records[order(names(records))]

# ── Build one cohort block: cohort-name strip on top, [violin | corr] below ──
make_cohort_block = function(rec, tag_violin, tag_corr) {
  v = rec$violin +
    labs(tag = tag_violin) +
    theme(plot.tag = element_text(face = "bold", size = 12),
          plot.tag.position = c(tag_x, tag_y))
  corr_p = rec$corr +
    labs(tag = tag_corr) +
    theme(plot.tag = element_text(face = "bold", size = 12),
          plot.tag.position = c(tag_x, tag_y))
  
  panel_row = wrap_plots(v, corr_p, ncol = 2, widths = c(w_violin, w_corr))
  
  title_strip = wrap_elements(full =
                                grid::textGrob(rec$name, gp = grid::gpar(fontface = "bold", fontsize = 16)))
  
  wrap_plots(title_strip, panel_row, ncol = 1,
             heights = c(title_strip_h, 1 - title_strip_h))
}

# ── Assemble: 2 cohort blocks per row, tags a..l ─────────────────────────────
tag_letters = letters[seq_len(2 * length(records))]
blocks = list()
for (i in seq_along(records)) {
  blocks[[length(blocks) + 1]] =
    make_cohort_block(records[[i]], tag_letters[2 * i - 1], tag_letters[2 * i])
}

combined = wrap_plots(blocks, ncol = 2)

n_cohorts = length(records)
n_rows    = ceiling(n_cohorts / 2)

ggsave("combined_violin_correlation_by_cohort.png",
       plot = combined,
       width = fig_width_in, height = row_height_in * n_rows,
       dpi = out_dpi, limitsize = FALSE)

# ════════════════════════════════════════════════════════════════════════════
# DCS target panel size — DETERMINISTIC.
#
# A row holds 2 cohort blocks; each block splits its width as w_violin : w_corr.
# So across the row the total relative width is 2 * (w_violin + w_corr); the
# violin panel gets w_violin of that, the corr panel w_corr.
#
# Panel HEIGHT = the block's panel-row height (row minus the title strip),
# matching what the live panels occupy. dcs_height_calib is the single manual
# nudge if it's still slightly off.
# ════════════════════════════════════════════════════════════════════════════
total_rel  = 2 * (w_violin + w_corr)
violin_w_in = fig_width_in * (w_violin / total_rel)
corr_w_in   = fig_width_in * (w_corr   / total_rel)

panel_row_h_in = row_height_in * (1 - title_strip_h) * dcs_height_calib

dcs_target = list(
  violin = list(width_in = violin_w_in, height_in = panel_row_h_in,
                panel_ar = ar_violin),
  corr   = list(width_in = corr_w_in, height_in = panel_row_h_in,
                panel_ar = ar_corr),
  dpi    = out_dpi
)
saveRDS(dcs_target, dcs_target_file)

message(sprintf(
  "Saved combined figure: %d cohorts, %d rows.\n  DCS violin PNG target: %.2f x %.2f in (panel_ar %.2f)\n  DCS corr   PNG target: %.2f x %.2f in (panel_ar %.2f)\n  Target -> %s. Copy to RIVM server, re-run the DCS script.",
  n_cohorts, n_rows,
  violin_w_in, panel_row_h_in, ar_violin,
  corr_w_in,   panel_row_h_in, ar_corr,
  dcs_target_file))
