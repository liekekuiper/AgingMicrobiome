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
from qiime2.plugins.feature_table.methods import filter_features_conditionally

# Select non-missing cases
def find_complete(metadata, model, subset, out, factors):
    # Ensure 'out' is a list
    if isinstance(out, str):
        out = [out]

    # Split the model into columns
    columns = model.split('+') if model else []

    # Combine columns and out
    all_columns = columns + out
    
    # Add specific columns if out is 'mortality'
    if out == ['mortality']:
        all_columns += ['studytime', 'age']

    # Drop rows with missing values in the specified columns
    meta_df = metadata.dropna(subset=all_columns)

    # Filter the DataFrame based on the subset
    if subset == 'all':
        final_df = meta_df
    elif subset in ['men', 'women']:
        final_df = meta_df[meta_df['sex'] == subset]
    elif subset.startswith('age_'):
        age_ranges = {
            'age_1': (18, 40),
            'age_2': (40, 50),
            'age_3': (50, 60),
            'age_4': (60, 70),
            'age_5': (70, float('inf'))
        }
        age_min, age_max = age_ranges[subset]
        final_df = meta_df[(meta_df['age'] >= age_min) & (meta_df['age'] < age_max)]
    else:
        raise ValueError(f"Invalid subset: {subset}")

    # Make a copy of the final_df to avoid SettingWithCopyWarning
    final_df = final_df.copy()

    # Determine factor names
    base_factors = ['ppump', 'metfor', 'statin', 'race']
    if factors != ['NA']:
        factor_names = base_factors + factors
    else:
        factor_names = base_factors

    if subset not in ['women', 'men']:
        factor_names = ['sex'] + factor_names

    # Re-level factors and convert to categorical
    for factor_name in factor_names:
        if factor_name == 'sex':
            final_df['sex'] = pd.Categorical(final_df['sex'], categories=['men', 'women'], ordered=True)
        elif factor_name == 'race':
            all_categories = final_df['race'].unique().tolist()
            if 'white' in all_categories:
                all_categories.remove('white')
            all_categories = ['white'] + all_categories
            final_df['race'] = pd.Categorical(final_df['race'], categories=all_categories, ordered=True)
        elif factor_name in factors:
            largest_category = final_df[factor_name].value_counts().idxmax()
            all_categories = [largest_category] + [cat for cat in final_df[factor_name].unique() if cat != largest_category]
            final_df[factor_name] = pd.Categorical(final_df[factor_name], categories=all_categories, ordered=True)
        else:
            final_df[factor_name] = pd.Categorical(final_df[factor_name], categories=[0, 1], ordered=True)

    # Create dummy variables for non-numeric categories if 'mortality' is in out
    if 'mortality' in out:
        for covariate in columns:
            if final_df[covariate].dtype == 'object' or final_df[covariate].dtype.name == 'category':
                final_df = pd.get_dummies(final_df, columns=[covariate], drop_first=True)

    return final_df

# Species-level: add unknown species when not known and replace spaces with _ for downstream analyses
def as_species(table, taxonomy):
    species = taxonomy['species'].to_dict()
    return table.collapse(
        lambda i, m: species.get(i, f'Unknown_Species_{i}').replace(' ', '_'), 
        norm=False, 
        axis='observation'
    )

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
        'shannon': ('alpha', table_ar, None, 'shannon_asv'),
        'chao1': ('alpha', table_ar, None, 'chao1_asv'),
        'simpson': ('alpha', table_ar, None, 'simpson_asv'),
        'simpson_e': ('alpha', table_ar, None, 'simpson_e_asv')
    }
    
    # Define the metrics and their respective parameters for the genus table
    genus_metrics = {
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

def process(taxonomy, tree, feature_table, output, threads, metadata, model, sub, out, factors):
    # Load metadata
    try:
        meta = pd.read_csv(metadata, sep='\t')
        meta.columns = [col.lower() for col in meta.columns]
        meta['sampleid'] = meta['sampleid'].astype(str) #Make strings for joining with feature table
        meta = meta.set_index('sampleid')
    except csv.Error as e:
        raise ValueError(f"Error parsing metadata file: {e}")
    meta_df = find_complete(meta, model, sub, out, factors)
    taxonomy = qiime2.Artifact.load(taxonomy).view(pd.DataFrame)
    tree_ar = qiime2.Artifact.load(tree)
    
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
    filtered_table_sample = feature_table.filter(common_sample_ids, axis='sample', inplace=False)
    filtered_sample_ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', filtered_table_sample)
    table_ar = filter_features_conditionally(filtered_sample_ar, abundance=0.01, prevalence=0.1).filtered_table

    #Genus level
    taxonomy['species'] = taxonomy['Taxon'].apply(lambda x: x.split('; ')[-1])
    species_table_tax = as_species(filtered_table_sample, taxonomy)
    species_table_ar_unfiltered = qiime2.Artifact.import_data('FeatureTable[Frequency]', species_table_tax)
    species_table_ar = filter_features_conditionally(species_table_ar_unfiltered, abundance=0.01, prevalence=0.1).filtered_table
    species_table = species_table_ar.view(biom.Table)
    species_table_clr = to_clr(species_table)

    # Calculate and add beta diversities to metadata
    meta_df = process_beta_diversities(table_ar, species_table_ar, tree_ar, threads, meta_df)
    
    # Calculate and add alpha diversities to metadata
    meta_df = process_alpha_diversities(table_ar, species_table_ar, tree_ar, threads, meta_df)
    speciestable_join = species_table_clr.transpose()
    newfile =  meta_df.join(speciestable_join)
    return newfile
if __name__ == '__main__':
    process()
