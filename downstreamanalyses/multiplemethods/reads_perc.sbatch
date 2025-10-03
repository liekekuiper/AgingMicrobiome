#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --time=1-00:00:00
#SBATCH -J genus_percentages
#SBATCH --mem=16G

# Load QIIME2 environment
source activate qiime2-2023.7

# Set cohort name here once!
COHORT="MrOS"

# Step 0: Create taxonomy table
qiime greengenes2 taxonomy-from-table \
  --i-reference-taxonomy ~/Microbiome/FirstTry/Deblur/2022.10.taxonomy.asv.nwk.qza \
  --i-table ${COHORT}.feature_table.qza \
  --o-classification ${COHORT}.taxonomy.qza

# Step 1: Group samples by 'preparation'
qiime feature-table group \
  --i-table ${COHORT}.feature_table.qza \
  --p-axis sample \
  --p-mode sum \
  --m-metadata-file ${COHORT}.metadata_with_sample.tsv \
  --m-metadata-column preparation \
  --o-grouped-table ${COHORT}.grouped_by_prep.qza

# Step 2: Collapse grouped table to genus level
qiime taxa collapse \
  --i-table ${COHORT}.grouped_by_prep.qza \
  --i-taxonomy ${COHORT}.taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ${COHORT}.genus_table.qza

# Step 3: Export collapsed table
EXPORT_DIR=genus_export_${COHORT}
mkdir -p ${EXPORT_DIR}
qiime tools export \
  --input-path ${COHORT}.genus_table.qza \
  --output-path ${EXPORT_DIR}

# Step 4: Convert BIOM to TSV
biom convert \
  -i ${EXPORT_DIR}/feature-table.biom \
  -o ${EXPORT_DIR}/${COHORT}.genus_table.tsv \
  --to-tsv

# Step 5: Calculate % reads per genus
python3 <<EOF
import pandas as pd

df = pd.read_csv("${EXPORT_DIR}/${COHORT}.genus_table.tsv", sep="\t", skiprows=1, index_col=0)
df = df.loc[:, ~df.columns.str.contains('^#')]  # Remove any metadata column

df_percent = df.div(df.sum(axis=0), axis=1) * 100
df_percent.to_csv("${EXPORT_DIR}/${COHORT}.genus_table_percentages.tsv", sep="\t")

print("Percentages saved to '${EXPORT_DIR}/${COHORT}.genus_table_percentages.tsv'")
EOF

