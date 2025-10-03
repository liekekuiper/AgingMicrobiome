import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import qiime2
import biom

### Settings
cohort_name = "MrOS"   # change per cohort
label_name  = "MrOS"   # label for figures and legends

### Input files
prev_file_16S = f"{cohort_name}_16S_Prev_Abundance_per_genus.csv"
prev_file_WGS = f"{cohort_name}_WGS_Prev_Abundance_per_genus.csv"
qza_file_16S  = f"{cohort_name}.16S.genus.qza"
qza_file_WGS  = f"{cohort_name}.WGS.genus.qza"

### Outputs
out_genus_csv   = f"{cohort_name}_per_genus_rank_correlations.csv"
out_sample_csv  = f"{cohort_name}_per_participant_rank_correlations.csv"
out_genus_plot  = f"{cohort_name}_per_genus_rank_correlation_density.png"
out_sample_plot = f"{cohort_name}_per_participant_rank_correlation_density.png"
legend_path     = f"Supp_Fig2_legend_{cohort_name}.txt"

### Helper functions
def summary_stats(series):
    """Return median and IQR (Q1, Q3)."""
    s = pd.Series(series).dropna()
    if s.empty:
        return np.nan, (np.nan, np.nan)
    med = s.median()
    q1 = s.quantile(0.25)
    q3 = s.quantile(0.75)
    return med, (q1, q3)

def save_density(df, col, main_title, outfile, n_participants=None, n_genera=None):
    """Create KDE density plot with titles showing N, genera, median ρ and IQR."""
    if df is None or df.empty:
        print(f"Skipping plot (no data): {outfile}")
        return

    med, (q1, q3) = summary_stats(df[col])

    parts = []
    if n_participants is not None:
        parts.append(f"N={n_participants}")
    if n_genera is not None:
        parts.append(f"Genera={n_genera}")
    parts.append(f"median ρ={med:.2f}")
    parts.append(f"IQR {q1:.2f}–{q3:.2f}")
    subtitle = ", ".join(parts)

    plt.figure(figsize=(9.2, 4.8))
    sns.kdeplot(df[col], fill=True)
    plt.xlabel("Spearman correlation (ρ)")
    plt.ylabel("Density")
    plt.xlim(0, 1)
    plt.axvline(med, linestyle="--")

    plt.suptitle(main_title, fontsize=12, y=1.02)
    plt.title(subtitle, fontsize=10)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved plot --> {outfile}")

### Prevalence filtering (≥10% at ≥1% relative abundance)
prev_16S = pd.read_csv(prev_file_16S)
prev_WGS = pd.read_csv(prev_file_WGS)

prev_16S_f = prev_16S[prev_16S["Percentage_with_≥0.01"] >= 10]
prev_WGS_f = prev_WGS[prev_WGS["Percentage_with_≥0.01"] >= 10]

common_taxa = set(prev_16S_f["Original_Taxon"]).intersection(prev_WGS_f["Original_Taxon"])
genera = sorted(common_taxa)
print(f"Shared genera passing prevalence threshold: {len(genera)}")
if len(genera) == 0:
    raise ValueError("No shared genera passing the prevalence threshold; cannot proceed.")

### Load genus-level tables from QIIME 2 artifacts
table_16S = qiime2.Artifact.load(qza_file_16S).view(biom.Table).to_dataframe()
table_WGS = qiime2.Artifact.load(qza_file_WGS).view(biom.Table).to_dataframe()

if hasattr(table_16S, "sparse"):
    table_16S = table_16S.sparse.to_dense()
if hasattr(table_WGS, "sparse"):
    table_WGS = table_WGS.sparse.to_dense()

### Subset to shared genera (rows = taxa)
table_16S = table_16S.loc[table_16S.index.intersection(genera)]
table_WGS = table_WGS.loc[table_WGS.index.intersection(genera)]

### Standardize sample IDs before transpose
table_16S.columns = table_16S.columns.str.replace(".16S", "", regex=False)
table_WGS.columns = table_WGS.columns.str.replace(".WGS", "", regex=False)

### Transpose to samples x genera
table_16S = table_16S.T
table_WGS = table_WGS.T

### Align samples and genera
common_samples = table_16S.index.intersection(table_WGS.index)
if len(common_samples) == 0:
    raise ValueError("No shared samples between 16S and WGS after ID harmonization.")
table_16S = table_16S.loc[common_samples, genera].copy()
table_WGS = table_WGS.loc[common_samples, genera].copy()

non_allzero = ~((table_16S.sum(axis=0) == 0) & (table_WGS.sum(axis=0) == 0))
table_16S = table_16S.loc[:, non_allzero]
table_WGS = table_WGS.loc[:, non_allzero]
genera = list(table_16S.columns)
print(f"Genera retained after zero-sum check: {len(genera)}")

### Per-genus correlations (across participants)
genus_results = []
for g in genera:
    x = table_16S[g]
    y = table_WGS[g]
    if x.var() == 0 or y.var() == 0:
        continue
    ρ, _ = spearmanr(x, y, nan_policy="omit")
    if np.isfinite(ρ):
        genus_results.append({"Genus": g, "Spearman_ρ": ρ})

genus_df = pd.DataFrame(genus_results)
if not genus_df.empty:
    genus_df.to_csv(out_genus_csv, index=False)
    print(f"Saved per-genus correlations: {out_genus_csv}")
else:
    print("No valid genera for correlation.")

### Per-participant correlations (across genera)
sample_results = []
X = table_16S[genera]
Y = table_WGS[genera]
for s in X.index:
    xv = X.loc[s]
    yv = Y.loc[s]
    if xv.var() == 0 or yv.var() == 0:
        continue
    ρ, _ = spearmanr(xv, yv, nan_policy="omit")
    if np.isfinite(ρ):
        sample_results.append({"Sample": s, "Spearman_ρ": ρ})

sample_df = pd.DataFrame(sample_results)
if not sample_df.empty:
    sample_df.to_csv(out_sample_csv, index=False)
    print(f"Saved per-participant correlations: {out_sample_csv}")
else:
    print("No valid samples for correlation.")

### Counts and summary stats
n_genera_input  = len(genera)
n_samples_input = len(common_samples)

n_genera_used   = len(genus_df)
n_samples_used  = len(sample_df)

g_med, (g_q1, g_q3) = summary_stats(genus_df["Spearman_ρ"]) if not genus_df.empty else (np.nan, (np.nan, np.nan))
s_med, (s_q1, s_q3) = summary_stats(sample_df["Spearman_ρ"]) if not sample_df.empty else (np.nan, (np.nan, np.nan))

### Save plots
save_density(
    genus_df, "Spearman_ρ",
    main_title=f"Per-genus Spearman correlations (16S vs metagenomics) — {label_name}",
    outfile=out_genus_plot,
    n_participants=n_samples_used,
    n_genera=n_genera_input
)

save_density(
    sample_df, "Spearman_ρ",
    main_title=f"Per-participant Spearman correlations (16S vs metagenomics) — {label_name}",
    outfile=out_sample_plot,
    n_participants=n_samples_used,
    n_genera=n_genera_input
)

### Write legend text file
legend_text = f"""Supplementary Fig. 2. Rank-based concordance between 16S rRNA gene amplicon and shotgun metagenomic sequencing in {label_name}.
(A) Per-genus Spearman correlations across participants (N participants used = {n_samples_used} of {n_samples_input} paired; Genera = {n_genera_input}); median ρ = {g_med:.2f}, IQR {g_q1:.2f}–{g_q3:.2f}.
(B) Per-participant Spearman correlations across genera (N participants used = {n_samples_used} of {n_samples_input} paired; Genera = {n_genera_input}); median ρ = {s_med:.2f}, IQR {s_q1:.2f}–{s_q3:.2f}.
Analyses were restricted to genera detected by both methods with prevalence ≥10% (relative abundance ≥1% in ≥10% of samples). Dashed lines indicate medians; density curves show the distributions of correlation coefficients."""
with open(legend_path, "w", encoding="utf-8") as f:
    f.write(legend_text)
print(f"Wrote legend --> {legend_path}")
