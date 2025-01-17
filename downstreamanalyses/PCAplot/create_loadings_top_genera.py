import click
import qiime2
import pandas as pd
import numpy as np
from skbio.stats.composition import clr
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import biom

# Function to aggregate data at the genus level using taxonomy
def as_genus(table, taxonomy):
    taxonomy['genus'] = taxonomy['Taxon'].apply(lambda x: x.split('; ')[-2])
    genus_map = taxonomy['genus'].to_dict()
    return table.collapse(lambda i, m: genus_map.get(i, 'Unknown'), norm=False, axis='observation')

# Function to calculate centered log-ratio (CLR) transformation
def to_clr(df):
    df += 1  # Add pseudocount
    df = df.div(df.sum(axis=0), axis=1)  # Convert to relative abundance
    return pd.DataFrame(clr(df.T), columns=df.index, index=df.columns).T

# Function to perform PCA
def perform_pca(data, n_components=3):
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(data)
    return pca, scores

@click.command()
@click.option('--taxonomy', type=click.Path(exists=True), required=True, help="Path to the taxonomy file.")
@click.option('--tables', type=click.Path(exists=True), multiple=True, required=True, help="Paths to qza feature tables for cohorts.")
@click.option('--cohort-names', multiple=True, required=True, help="Names of cohorts corresponding to the tables.")
@click.option('--output-dir', type=click.Path(exists=False), required=True, help="Directory to save outputs.")
@click.option('--top-n', type=int, default=5, help="Number of top genera to filter.")
def get_loadings(taxonomy, tables, cohort_names, output_dir, top_n):
    if len(tables) != len(cohort_names):
        raise ValueError("Number of tables must match the number of cohort names.")
    
    # Load taxonomy
    taxonomy_df = qiime2.Artifact.load(taxonomy).view(pd.DataFrame)
    
    # Load and process BIOM tables
    biom_tables = [qiime2.Artifact.load(tbl).view(biom.Table) for tbl in tables]
    genus_tables = [as_genus(tbl, taxonomy_df) for tbl in biom_tables]

    # Identify common top genera
    top_genera_per_cohort = [
        set(pd.Series(table.sum(axis="observation"), index=table.ids(axis="observation"))
            .sort_values(ascending=False)
            .head(top_n)
            .index.tolist())
        for table in genus_tables]
    common_top_genera = set.intersection(*top_genera_per_cohort)
    
    # Save common_top_genera to a text file
    with open(f"{output_dir}/common_genera.txt", "w") as f:
        for genus in common_top_genera:
            f.write(f"{genus}\n")


    # Prepare data for combined PCA
    cohort_labels = []
    combined_clr_data = []

    for table, cohort_name in zip(genus_tables, cohort_names):
        genus_table_filtered = table.filter(common_top_genera, axis='observation').to_dataframe()
        clr_data = to_clr(genus_table_filtered.T)
        clr_data["Cohort"] = cohort_name
        cohort_labels.extend([cohort_name] * clr_data.shape[0])
        combined_clr_data.append(clr_data.drop(columns=["Cohort"]))

    # Combine all CLR-transformed data
    combined_clr_data = pd.concat(combined_clr_data, axis=0)
    cohort_labels = pd.Series(cohort_labels, index=combined_clr_data.index, name="Cohort")

    # Perform PCA on combined data
    pca, pca_scores = perform_pca(combined_clr_data, n_components=3)
    loadings = pd.DataFrame(pca.components_.T, index=combined_clr_data.columns, columns=[f"PC{j+1}" for j in range(3)])
    scores = pd.DataFrame(pca_scores, index=combined_clr_data.index, columns=[f"PC{j+1}" for j in range(3)])
    scores["Cohort"] = cohort_labels  # Add cohort labels

    # Save PCA results
    scores.to_csv(f"{output_dir}/combined_pca_scores.csv", index=True)
    loadings.to_csv(f"{output_dir}/combined_pca_loadings.csv", index=True)

    # Print explained variance
    explained_variance = pca.explained_variance_ratio_
    # After PCA computation
    explained_variance_series = pd.Series(explained_variance, index=[f"PC{i+1}" for i in range(len(explained_variance))])
    explained_variance_series.to_csv(f"{output_dir}/explained_variance.csv")

    print("Explained variance by each principal component:", explained_variance)

if __name__ == '__main__':
    get_loadings()
