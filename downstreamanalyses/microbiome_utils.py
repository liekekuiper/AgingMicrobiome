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
    if isinstance(out, str):
        out = [out]

    columns = model.split('+') if model else []
    all_columns = columns + out

    if out == ['mortality']:
        all_columns += ['studytime', 'age']

    meta_df = metadata.dropna(subset=all_columns)

    # Restrict to participants aged 18 and older
    meta_df = meta_df[meta_df['age'] >= 18]

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

    final_df = final_df.copy()
    base_factors = ['ppump', 'metfor', 'statin', 'race']
    factor_names = base_factors + factors if factors != ['NA'] else base_factors
    if subset not in ['women', 'men']:
        factor_names = ['sex'] + factor_names

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

    if 'mortality' in out:
        for covariate in columns:
            if final_df[covariate].dtype == 'object' or final_df[covariate].dtype.name == 'category':
                final_df = pd.get_dummies(final_df, columns=[covariate], drop_first=True)

    return final_df

def as_genus(table, taxonomy):
    genus = taxonomy['genus'].to_dict()
    return table.collapse(lambda i, m: genus.get(i, f'Unknown_Genus_{i}'), norm=False, axis='observation')

def to_clr(data):
    df = data.to_dataframe()
    df += 1
    df = df.div(df.sum(axis=0), axis=1)
    return pd.DataFrame(clr(df.T), columns=df.index, index=df.columns).T

def calculate_min_dissimilarity(distance_matrix):
    temp_df = distance_matrix.view(DistanceMatrix)
    dm_df = temp_df.to_data_frame()
    dm_matrix = dm_df.values
    np.fill_diagonal(dm_matrix, np.nan)
    return np.nanmin(dm_matrix, axis=1)

def add_alpha_diversity_to_metadata(metadata_df, diversity_metric, column_name):
    alpha_df = diversity_metric.view(pd.Series)
    metadata_df[column_name] = metadata_df.index.map(alpha_df)
    return metadata_df

def process_beta_diversities(table_ar, genus_table_ar, species_table_ar, tree_ar, threads, metadata):
    beta_metrics = {
        'braycurtis': ['min_bray_asv', 'min_bray_genus'],
        'jaccard': ['min_jacc_asv', 'min_jacc_genus']
    }

    # Add species-level metrics if applicable
    if species_table_ar is not None:
        beta_metrics['braycurtis'].append('min_bray_species')
        beta_metrics['jaccard'].append('min_jacc_species')

    for metric, columns in beta_metrics.items():
        print(metric)

        # ASV-level (only if applicable, i.e., not None)
        if columns[0] is not None:
            dm_asv = diversity.actions.beta(table_ar, metric=metric).distance_matrix
            metadata[columns[0]] = calculate_min_dissimilarity(dm_asv)

        # Genus-level
        dm_genus = diversity.actions.beta(genus_table_ar, metric=metric).distance_matrix
        metadata[columns[1]] = calculate_min_dissimilarity(dm_genus)

        # Species-level (if defined)
        if species_table_ar is not None and len(columns) > 2:
            dm_species = diversity.actions.beta(species_table_ar, metric=metric).distance_matrix
            metadata[columns[2]] = calculate_min_dissimilarity(dm_species)


    print('uu')
    dm_uu = diversity.actions.beta_phylogenetic(table_ar, tree_ar, threads=threads, metric='unweighted_unifrac').distance_matrix
    metadata['min_uu_feature'] = calculate_min_dissimilarity(dm_uu)

    print('wu')
    dm_wu = diversity.actions.beta_phylogenetic(table_ar, tree_ar, threads=threads, metric='weighted_normalized_unifrac').distance_matrix
    metadata['min_wu_feature'] = calculate_min_dissimilarity(dm_wu)

    return metadata

def process_alpha_diversities(table_ar, genus_table_ar, species_table_ar, tree_ar, threads, metadata):
    all_metrics = {
        'asv': (table_ar, 'asv'),
        'genus': (genus_table_ar, 'genus')
    }

    if species_table_ar is not None:
        all_metrics['species'] = (species_table_ar, 'species')

    metric_names = ['shannon', 'chao1', 'simpson', 'simpson_e']

    def process_metric(method, table, tree, metric, column_name, metadata):
        if method == 'alpha_phylogenetic':
            dm = getattr(diversity.actions, method)(table, tree, metric=metric)
        else:
            dm = getattr(diversity.actions, method)(table, metric=metric)
        return add_alpha_diversity_to_metadata(metadata, dm.alpha_diversity, column_name)

    for table_label, (table, suffix) in all_metrics.items():
        for metric in metric_names:
            method = 'alpha'
            tree = None
            column_name = f'{metric}_{suffix}'
            print((method, table, tree, column_name))
            metadata = process_metric(method, table, tree, metric, column_name, metadata)

    return metadata

def process(taxonomy, tree, feature_table, output, threads, metadata, model, sub, out, factors, label='16s'):
    try:
        meta = pd.read_csv(metadata, sep='\t')
        meta.columns = [col.lower() for col in meta.columns]
        meta['sampleid'] = meta['sampleid'].astype(str)
        meta = meta.set_index('sampleid')
    except csv.Error as e:
        raise ValueError(f"Error parsing metadata file: {e}")

    meta_df = find_complete(meta, model, sub, out, factors)
    taxonomy = qiime2.Artifact.load(taxonomy).view(pd.DataFrame)
    tree_ar = qiime2.Artifact.load(tree)
    sample_ids = meta_df.index.tolist()

    ftable_ar = qiime2.Artifact.load(feature_table)
    feature_table = ftable_ar.view(biom.Table)
    valid_sample_ids = set(feature_table.ids(axis='sample'))
    common_sample_ids = [sid for sid in sample_ids if sid in valid_sample_ids]
    meta_df = meta_df.loc[common_sample_ids]
    filtered_table_sample = feature_table.filter(common_sample_ids, axis='sample', inplace=False)
    filtered_sample_ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', filtered_table_sample)

    taxonomy['genus'] = taxonomy['Taxon'].apply(lambda x: x.split('; ')[-2])
    genus_table_tax = as_genus(filtered_table_sample, taxonomy)
    genus_table_ar_unfiltered = qiime2.Artifact.import_data('FeatureTable[Frequency]', genus_table_tax)
    genus_table_ar = filter_features_conditionally(genus_table_ar_unfiltered, abundance=0.01, prevalence=0.1).filtered_table
    genus_table = genus_table_ar.view(biom.Table)
    genus_table_clr = to_clr(genus_table)

    species_table_ar_unfiltered = None
    if label.lower() != '16s':
        print('Calculating species-level metrics...')
        taxonomy['species'] = taxonomy['Taxon'].apply(lambda x: x.split('; ')[-1])
        species_table_tax = filtered_table_sample.collapse(
            lambda i, m: taxonomy['species'].get(i, f'Unknown_Species_{i}'),
            norm=False,
            axis='observation'
        )
        species_table_ar_unfiltered = qiime2.Artifact.import_data('FeatureTable[Frequency]', species_table_tax)
        species_table_ar = filter_features_conditionally(species_table_ar_unfiltered, abundance=0.01, prevalence=0.1).filtered_table

    meta_df = process_beta_diversities(
        filtered_sample_ar,
        genus_table_ar_unfiltered,
        species_table_ar_unfiltered,
        tree_ar,
        threads,
        meta_df
    )

    meta_df = process_alpha_diversities(
        filtered_sample_ar,
        genus_table_ar_unfiltered,
        species_table_ar_unfiltered,
        tree_ar,
        threads,
        meta_df
    )

    genustable_join = genus_table_clr.transpose()
    newfile = meta_df.join(genustable_join)

    if species_table_ar is not None:
        species_table = species_table_ar.view(biom.Table)
        species_table_clr = to_clr(species_table)
        species_table_join = species_table_clr.transpose()
        newfile = newfile.join(species_table_join)

    return newfile


if __name__ == '__main__':
    process()
