import pandas as pd
import numpy as np
import skbio
import biom
import qiime2
from qiime2.plugins import diversity


### Change cohort name
cohort = 'MrOS' #Change cohort name

threads=6 #number of threads
tree="2022.10.phylogeny.asv.nwk.qza"
tree_ar = qiime2.Artifact.load(tree)


#### Define functions
# Genus-level
def as_genus(table, taxonomy):
    genus = taxonomy['genus'].to_dict()
    return table.collapse(lambda i, m: genus.get(i, f'Unknown_Genus_{i}'), norm=False, axis='observation')

# Function to add alpha diversity to metadata
def add_alpha_diversity_to_metadata(metadata_df, diversity_metric, column_name):
    alpha_df = diversity_metric.view(pd.Series)
    metadata_df[column_name] = metadata_df.index.map(alpha_df)
    return metadata_df

#Function to calculate alpha diversity
def process_alpha_diversities(table_ar, genus_table_ar, tree_ar, threads, metadata):
    # Combined metrics with table and column suffix
    all_metrics = {
        'asv': (table_ar, 'asv'),
        'genus': (genus_table_ar, 'genus')
    }

    metric_names = ['shannon', 'chao1', 'simpson', 'simpson_e']

    # Function to process a single metric
    def process_metric(method, table, tree, metric, column_name, metadata):
        if method == 'alpha_phylogenetic':
            dm = getattr(diversity.actions, method)(table, tree, metric=metric)
        else:
            dm = getattr(diversity.actions, method)(table, metric=metric)
        return add_alpha_diversity_to_metadata(metadata, dm.alpha_diversity, column_name)

    # Loop through each table and its suffix
    for table_label, (table, suffix) in all_metrics.items():
        for metric in metric_names:
            method = 'alpha'
            tree = None
            column_name = f'{metric}_{suffix}'
            print((method, table, tree, column_name))
            metadata = process_metric(method, table, tree, metric, column_name, metadata)

    return metadata

feature_table = qiime2.Artifact.load(f'{cohort}.feature_table.qza')
taxonomy = qiime2.Artifact.load(f'{cohort}.taxonomy.qza').view(pd.DataFrame)
taxonomy['genus'] = taxonomy['Taxon'].apply(lambda x: x.split('; ')[-2])

filtered_table_sample = feature_table.view(biom.Table)
genus_table_tax = as_genus(filtered_table_sample, taxonomy)
genus_table_ar_unfiltered = qiime2.Artifact.import_data('FeatureTable[Frequency]', genus_table_tax)

meta = pd.read_csv(f'{cohort}.metadata_with_sample.tsv', sep='\t')
meta.columns = [col.lower() for col in meta.columns]
meta['#sampleid'] = meta['#sampleid'].astype(str) #Make strings for joining with feature table
meta = meta.set_index('#sampleid')
meta.index.names = ['sampleid']


meta_df = process_alpha_diversities(feature_table, genus_table_ar_unfiltered, tree_ar, threads, meta)
# Select only relevant columns
alpha_cols = ['shannon_asv', 'chao1_asv', 'simpson_asv', 'simpson_e_asv',
              'shannon_genus', 'chao1_genus', 'simpson_genus', 'simpson_e_genus']

# Pivoting to wide format based on sample and preparation
wide_df = meta_df.pivot(index='sample', columns='preparation', values=alpha_cols)

# Flatten multi-level columns
wide_df.columns = [f'{col}_{prep}' for col, prep in wide_df.columns]

# Get list of unique preparations
preparations = meta_df['preparation'].unique()

# Create list to store correlation results
correlation_results = []

# Loop through all alpha metrics
for metric in alpha_cols:
    # Get columns for this metric across preparations
    metric_cols = [f"{metric}_{prep}" for prep in preparations if f"{metric}_{prep}" in wide_df.columns]
    
    # Compute pairwise correlations
    for i in range(len(metric_cols)):
        for j in range(i + 1, len(metric_cols)):
            col1 = metric_cols[i]
            col2 = metric_cols[j]
            corr = wide_df[[col1, col2]].corr(method = 'spearman').iloc[0, 1]
            correlation_results.append({
                'metric': metric,
                'prep1': col1.split('_')[-1],
                'prep2': col2.split('_')[-1],
                'correlation': corr
            })

# Convert to DataFrame
cor_df = pd.DataFrame(correlation_results)


cor_df.to_csv('correlation_alpha16SWGS.csv')
