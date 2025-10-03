#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -J abundance
#SBATCH --time=1-00:00:00

source activate qiime2-2023.7
set -e
set -x

# ---- CONFIG ----
COHORT="MrOS"  # <<< Only change this!
TAXONOMY_REF="~/Microbiome/FirstTry/Deblur/2022.10.taxonomy.asv.nwk.qza"
META="${COHORT}.metadata_with_sample.tsv"
TABLE="${COHORT}.feature_table.qza"

# ---- Loop over preparation types ----
for PREP in "16S" "WGS"; do

  # Step 1: Filter feature table by preparation type
  qiime feature-table filter-samples \
    --i-table ${TABLE} \
    --m-metadata-file ${META} \
    --p-where "[preparation] = '${PREP}'" \
    --o-filtered-table ${COHORT}.${PREP}.filtered.qza

  # Step 2: Generate taxonomy for this filtered table
  qiime greengenes2 taxonomy-from-table \
    --i-reference-taxonomy ~/Microbiome/FirstTry/Deblur/2022.10.taxonomy.asv.nwk.qza \
    --i-table ${COHORT}.${PREP}.filtered.qza \
    --o-classification ${COHORT}.${PREP}.taxonomy.qza

  # Step 3: Collapse to genus level
  qiime taxa collapse \
    --i-table ${COHORT}.${PREP}.filtered.qza \
    --i-taxonomy ${COHORT}.${PREP}.taxonomy.qza \
    --p-level 6 \
    --o-collapsed-table ${COHORT}.${PREP}.genus.qza

  # Step 4: Convert to relative abundance
  qiime feature-table relative-frequency \
    --i-table ${COHORT}.${PREP}.genus.qza \
    --o-relative-frequency-table ${COHORT}.${PREP}.rel_genus.qza

  # Step 5: Python to calculate prevalence/abundance
  python3 <<EOF
import pandas as pd
from qiime2 import Artifact
import biom

rel_genus = Artifact.load("${COHORT}.${PREP}.rel_genus.qza").view(biom.Table).to_dataframe()

def extract_genus(taxon):
    parts = taxon.split(';')
    for part in reversed(parts):
        if part.startswith('g__'):
            return part
    return taxon

genus_names = rel_genus.index.to_series().map(extract_genus)

percent_with_001 = (rel_genus >= 0.01).sum(axis=1) / rel_genus.shape[1] * 100
mean_abundance = rel_genus.mean(axis=1)
median_abundance = rel_genus.median(axis=1)
sd_abundance = rel_genus.std(axis=1)

prev_df = pd.DataFrame({
    'Original_Taxon': rel_genus.index,
    'Genus': genus_names,
    'Percentage_with_≥0.01': percent_with_001,
    'Mean_Abundance': mean_abundance,
    'Median_Abundance': median_abundance,
    'SD_Abundance': sd_abundance
})

prev_df.to_csv("${COHORT}_${PREP}_Prev_Abundance_per_genus.csv", index=False)
print(f"Saved: ${COHORT}_${PREP}_Prev_Abundance_per_genus.csv")
EOF

done

