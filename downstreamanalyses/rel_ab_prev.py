import pandas as pd
import numpy as np
from microbiome_utils import *

# Load relative table
rel_genus = qiime2.Artifact.load('rel-genus-table.qza').view(biom.Table).to_dataframe()

# Get genus name
def extract_genus(taxon):
    parts = taxon.split(';')
    for part in reversed(parts):
        if part.startswith('g__'):
            return part
    return taxon

# Map genus name
genus_names = rel_genus.index.to_series().map(extract_genus)

# Calculate percentage participants abundance ≥ 0.01
percent_with_001 = (rel_genus >= 0.01).sum(axis=1) / rel_genus.shape[1] * 100

mean_abundance = rel_genus.mean(axis=1)
median_abundance = rel_genus.median(axis=1)
sd_abundance = rel_genus.std(axis=1)

# Save results
prev_df = pd.DataFrame({
    'Original_Taxon': rel_genus.index,
    'Genus': genus_names,
    'Percentage_with_≥0.01': percent_with_001,
    'Mean_Abundance': mean_abundance,
    'Median_Abundance': median_abundance,
    'SD_Abundance': sd_abundance
})

prev_df.to_csv('Prev_Abundance_per_genus.csv', index = False)
