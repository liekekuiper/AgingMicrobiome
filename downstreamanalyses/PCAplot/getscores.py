import pandas as pd
import numpy as np
from skbio.stats.composition import clr
import qiime2
import biom

cohort_name = '' #Type cohort name if the cohort has both 16S and shotgun make either cohortA16S or cohortAshotgun
feature_table_path = '' #Type the path to the feature table (feature.table.gg2-2022.10.qza)
taxonomy_path = '' #Type the path to the taxonomy table (df.gg2.taxonomy.qza)
pca_loadings_path = 'combined_pca_loadings.csv' #Type the path to file with the pca loadings based on FHS (16S & shotgun), MrOS (16S & shotgun), and Rotterdam Study (16S) data
common_genera_path = 'common_genera.txt' #Type the path to the file with the 4 genera amongst the top 10 genera in FHS (16S & shotgun), MrOS (16S & shotgun), and Rotterdam Study (16S) data


####### Keep everything underneath this line untouched ########

# Function to calculate centered log-ratio (CLR) transformation
def to_clr(df):
    df += 1  # Add pseudocount
    df = df.div(df.sum(axis=0), axis=1)  # Convert to relative abundance
    return pd.DataFrame(clr(df.T), columns=df.index, index=df.columns).T

# Function to compute PCA scores using predefined loadings
def apply_loadings(new_data, loadings):
    return new_data.dot(loadings)

#On genus level
def as_genus(table, taxonomy):
    taxonomy['genus'] = taxonomy['Taxon'].apply(lambda x: x.split('; ')[-2])
    genus_map = taxonomy['genus'].to_dict()
    return table.collapse(lambda i, m: genus_map.get(i, 'Unknown'), norm=False, axis='observation')

# Load loadings for PCA
loadings = pd.read_csv(pca_loadings_path, index_col=0)

# Load taxonomy and featuretable QIIME 2 artifacts
taxonomy_artifact = qiime2.Artifact.load(taxonomy_path)
taxonomy_df = taxonomy_artifact.view(pd.DataFrame)

feature_table_artifact =qiime2.Artifact.load(feature_table_path).view(biom.Table)
feature_table = as_genus(feature_table_artifact, taxonomy_df)

# Filter the table for common genera
common_top_genera = pd.read_csv(common_genera_path, header=None).squeeze("columns").tolist()
genus_table_filtered = feature_table.filter(common_top_genera, axis='observation').to_dataframe()

# Perform CLR transformation
clr_data = to_clr(genus_table_filtered.T)

# Apply loadings to compute scores
scores = apply_loadings(clr_data, loadings)

#Make index anonymus for sharing
scores.index = [f"{i+1}_{cohort_name}" for i in range(scores.shape[0])]

# Add the cohort name as a new column
scores["Cohort"] = cohort_name

# Save the scores
scores.to_csv(f"pca_scores_{cohort_name}_common_genera.csv")

print(f"PCA scores for the new cohort '{cohort_name}' have been saved.")
