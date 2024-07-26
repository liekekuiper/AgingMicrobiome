import click
import qiime2
from qiime2.plugins import diversity
import pandas as pd
import numpy as np
from skbio.stats.composition import clr
import skbio
import biom
import csv
from skbio import DistanceMatrix


# Select non-missing cases
def find_complete(metadata, model, subset):
    # Extract the column names from the model string
    columns = model.split('+')
    # Drop rows with missing values in the specified columns
    meta_df = metadata.dropna(subset=columns)
    
    if subset == 'all':
        final_df = meta_df
    elif subset == 'men':
        final_df = meta_df[meta_df['sex'] == 'men']
    elif subset == 'women':
        final_df = meta_df[meta_df['sex'] == 'women']
    elif subset == 'age_1':
        final_df = meta_df[(meta_df['age'] >= 18) & (meta_df['age'] < 40)]
    elif subset == 'age_2':
        final_df = meta_df[(meta_df['age'] >= 40) & (meta_df['age'] < 50)]
    elif subset == 'age_3':
        final_df = meta_df[(meta_df['age'] >= 50) & (meta_df['age'] < 60)]
    elif subset == 'age_4':
        final_df = meta_df[(meta_df['age'] >= 60) & (meta_df['age'] < 70)]
    elif subset == 'age_5':
        final_df = meta_df[(meta_df['age'] >= 70)]

# Determine factor names
    if subset in ['women', 'men']:
        factor_names = ['ppump', 'metfor', 'statin', 'race']
    else:
        factor_names = ['sex', 'ppump', 'metfor', 'statin', 'race']
    
    # Re-level factors and convert to categorical
    for factor_name in factor_names:
        if factor_name == 'sex':
            final_df['sex'] = pd.Categorical(final_df['sex'], categories=['men', 'women'])
            final_df['sex'] = final_df['sex'].cat.reorder_categories(['men', 'women'], ordered=True)
        elif factor_name == 'race':
            # Ensure 'white' is the reference category
            all_categories = final_df['race'].unique().tolist()
            if 'white' in all_categories:
                all_categories.remove('white')
            all_categories = ['white'] + all_categories
            final_df['race'] = pd.Categorical(final_df['race'], categories=all_categories, ordered=True)
        else:
            final_df[factor_name] = pd.Categorical(final_df[factor_name], categories=['0', '1'])
            final_df[factor_name] = final_df[factor_name].cat.reorder_categories(['0', '1'], ordered=True)
    
    return final_df


# Genus-level
def as_genus(table, taxonomy):
    genus = taxonomy['genus'].to_dict()
    return table.collapse(lambda i, m: genus.get(i, f'Unknown Genus ({i})'), norm=False, axis='observation')


# Function from absolute abundances to clr
def to_clr(data):
    df = data.to_dataframe()
    df += 1                              # add pseudocount
    df = df.div(df.sum(axis=0), axis=1)  # relative abundance
    return pd.DataFrame(clr(df.T), columns=df.index, index=df.columns).T # clr

# Function to calculate minimum dissimilarity (uniqueness)
def calculate_min_dissimilarity(distance_matrix):
    temp_df = distance_matrix.view(DistanceMatrix)
    dm_df = temp_df.to_data_frame()
    dm_matrix = dm_df.values
    np.fill_diagonal(dm_matrix, np.nan)  # Replace diagonal with NaNs to ignore self-comparisons
    return np.nanmin(dm_matrix, axis=1)

# Function to add alpha diversity to metadata
def add_alpha_diversity_to_metadata(metadata_df, diversity_metric, column_name):
    alpha_df = diversity_metric.view(pd.Series)
    metadata_df[column_name] = metadata_df.index.map(alpha_df)
    return metadata_df

# Helper function to process beta diversities
def process_beta_diversities(table_ar, genus_table_ar, tree_ar, threads, metadata):
    beta_metrics = {
        'braycurtis': ['min_bray_asv', 'min_bray_genus'],
        'jaccard': ['min_jacc_asv', 'min_jacc_genus'],
        'aitchison': ['min_aitch_asv', 'min_aitch_genus'] #diversity.actions.beta does a clr in the function so we use non-clr transformed data
    }
    for metric, columns in beta_metrics.items():
        print(metric)
        dm_asv = diversity.actions.beta(table_ar, metric=metric).distance_matrix
        metadata[columns[0]] = calculate_min_dissimilarity(dm_asv)
        dm_genus = diversity.actions.beta(genus_table_ar, metric=metric).distance_matrix
        metadata[columns[1]] = calculate_min_dissimilarity(dm_genus)
    print('uu')
    dm_uu = diversity.actions.beta_phylogenetic(table_ar, tree_ar, threads=threads, metric='unweighted_unifrac').distance_matrix
    metadata['min_uu_feature'] = calculate_min_dissimilarity(dm_uu)
    print('wu')
    dm_wu = diversity.actions.beta_phylogenetic(table_ar, tree_ar, threads=threads, metric='weighted_normalized_unifrac').distance_matrix
    metadata['min_wu_feature'] = calculate_min_dissimilarity(dm_wu)
    
    return metadata

def process_alpha_diversities(table_ar, genus_table_ar, tree_ar, threads, metadata):
    # Define the metrics and their respective parameters for the main table
    alpha_metrics = {
#       'faith_pd': ('alpha_phylogenetic', table_ar, tree_ar, 'faith_pd_asv'),
        'shannon': ('alpha', table_ar, None, 'shannon_asv'),
        'chao1': ('alpha', table_ar, None, 'chao1_asv'),
        'simpson': ('alpha', table_ar, None, 'simpson_asv'),
        'simpson_e': ('alpha', table_ar, None, 'simpson_e_asv')
    }
    
    # Define the metrics and their respective parameters for the genus table
    genus_metrics = {
#        'faith_pd': ('alpha_phylogenetic', genus_table_ar, tree_ar, 'faith_pd_genus'),
        'shannon': ('alpha', genus_table_ar, None, 'shannon_genus'),
        'chao1': ('alpha', genus_table_ar, None, 'chao1_genus'),
        'simpson': ('alpha', genus_table_ar, None, 'simpson_genus'),
        'simpson_e': ('alpha', genus_table_ar, None, 'simpson_e_genus')
    }

    # Function to process a single metric
    def process_metric(method, table, tree, metric, column_name, metadata):
        if method == 'alpha_phylogenetic':
            dm = getattr(diversity.actions, method)(table, tree,  metric=metric)
        else:
            dm = getattr(diversity.actions, method)(table, metric=metric)
        return add_alpha_diversity_to_metadata(metadata, dm.alpha_diversity, column_name)

    # Process all alpha metrics for the asv table
    for metric, params in alpha_metrics.items():
        print(params)
        method, table, tree, column_name = params
        metadata = process_metric(method, table, tree, metric, column_name, metadata)

    # Process all alpha metrics for the genus table
    for metric, params in genus_metrics.items():
        print(params)
        method, table, tree, column_name = params
        metadata = process_metric(method, table, tree, metric, column_name, metadata)
    
    return metadata


#@click.command()
#@click.option('--taxonomy', type=click.Path(exists=True), required=True)
#@click.option('--tree', type=click.Path(exists=True), required=True)
#@click.option('--feature-table', type=click.Path(exists=True), required=True)
#@click.option('--metadata', type=click.Path(exists=True), required=True)
#@click.option('--output', type=click.Path(exists=False), required=True)
#@click.option('--threads', type=int, required=True)
def process(taxonomy, tree, feature_table, output, threads, metadata, model, sub):
    # Load metadata
    try:
        meta = pd.read_csv(metadata, sep='\t', dtype=str)
        meta.columns = map(str.lower, meta.columns)
        meta = meta.set_index('sampleid')
    except csv.Error as e:
        raise ValueError(f"Error parsing metadata file: {e}")
    meta_df = find_complete(meta, model, sub)
    taxonomy = qiime2.Artifact.load(taxonomy).view(pd.DataFrame)
    tree_ar = qiime2.Artifact.load(tree)
    tree = tree_ar.view(skbio.TreeNode)
    
    #Save sample ids
    sample_ids = meta_df.index.tolist()
    
    #Load data
    ftable_ar = qiime2.Artifact.load(feature_table)
    feature_table = ftable_ar.view(biom.Table)
    
    # Filter table to keep only sample IDs present in metadata and feature table
    # Get sample IDs from feature table
    valid_sample_ids = set(feature_table.ids(axis='sample'))

    # Filter sample IDs to keep only those present in both meta_df and table_16s
    common_sample_ids = [sid for sid in sample_ids if sid in valid_sample_ids]
    meta_df = meta_df.loc[common_sample_ids]
    filtered_table = feature_table.filter(common_sample_ids, axis='sample', inplace=False)

    #Make Qiime object again
    table_ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', filtered_table)
    #Perform CLR
    table_clr = to_clr(filtered_table)
    table_clr_ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', table_clr)

    #Genus level
    taxonomy['genus'] = taxonomy['Taxon'].apply(lambda x: x.split('; ')[-2])
    genus_table = as_genus(filtered_table, taxonomy)
    genus_table_ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', genus_table)
    genus_table_clr = to_clr(genus_table)
    genus_table_clr_ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', genus_table_clr)
#    genus_table_clr_ar.save(output + '.clr.genus.feature_table.qza')

    # Calculate and add beta diversities to metadata
    meta_df = process_beta_diversities(table_ar, genus_table_ar, tree_ar, threads, meta_df)
    
    # Calculate and add alpha diversities to metadata
    meta_df = process_alpha_diversities(table_ar, genus_table_ar, tree_ar, threads, meta_df)

    meta_df.to_csv(output + '.metadata.tsv', sep='\t', index=True, header=True)
    genustable_join = genus_table_clr.transpose()
    newfile =  meta_df.join(genustable_join)
    newfile.to_csv(output + '.checkfile.tsv', sep='\t', index=True, header=True)
    return newfile
if __name__ == '__main__':
    process()
