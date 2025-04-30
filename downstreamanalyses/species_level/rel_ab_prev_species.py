import pandas as pd
import numpy as np
from microbiome_utils import *

# Load relative table
rel_species = qiime2.Artifact.load('rel-species-table.qza').view(biom.Table).to_dataframe()

# Get species name
def extract_species(taxon):
    parts = taxon.split(';')
    for part in reversed(parts):
        if part.startswith('s__'):
            return part
    return taxon

# Map species name
species_names = rel_species.index.to_series().map(extract_species)

# Calculate percentage participants abundance ≥ 0.01
percent_with_001 = (rel_species >= 0.01).sum(axis=1) / rel_species.shape[1] * 100

mean_abundance = rel_species.mean(axis=1)
median_abundance = rel_species.median(axis=1)
sd_abundance = rel_species.std(axis=1)

# Save results
prev_df = pd.DataFrame({
    'Original_Taxon': rel_species.index,
    'Species': species_names,
    'Percentage_with_≥0.01': percent_with_001,
    'Mean_Abundance': mean_abundance,
    'Median_Abundance': median_abundance,
    'SD_Abundance': sd_abundance
})

prev_df.to_csv('Prev_Abundance_per_species.csv', index = False)
