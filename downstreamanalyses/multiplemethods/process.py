import click
import qiime2
from qiime2.plugins import diversity, emperor
import pandas as pd
import skbio
import biom


def as_genus(table, taxonomy):
    genus = taxonomy['genus'].to_dict()

    return table.collapse(lambda i, m: genus[i], norm=False,
                          axis='observation')

def concat(table_16s, table_wgs, metadata):
    md_ids = set(metadata.index)
    overlap = (set(table_16s.ids()) & set(table_wgs.ids()) & md_ids)
    table_16s = table_16s.filter(overlap).remove_empty()
    table_wgs = table_wgs.filter(overlap).remove_empty()
    table_16s.update_ids({i: i + '.16S' for i in table_16s.ids()}, inplace=True)
    table_wgs.update_ids({i: i + '.WGS' for i in table_wgs.ids()}, inplace=True)
    table = table_16s.concat(table_wgs)

    metadata = metadata.loc[overlap]
    md16s = metadata.copy()
    mdwgs = metadata.copy()
    md16s['preparation'] = '16S'
    mdwgs['preparation'] = 'WGS'
    md16s.index = [i + '.16S' for i in md16s.index]
    mdwgs.index = [i + '.WGS' for i in mdwgs.index]
    md16s.index.name = '#SampleID'
    mdwgs.index.name = '#SampleID'
    md = pd.concat([md16s, mdwgs])

    return table, md


#def subsample(table, tree, depth):
#    feats = {n.name for n in tree.tips()}
#    table = table.filter(feats & set(table.ids(axis='observation')), axis='observation')
#    table = table.filter(lambda v, i, m: v.sum() >= depth).remove_empty()
#    return table.subsample(depth, with_replacement=True)

def filter_table(table, tree, depth=None):
    feats = {n.name for n in tree.tips()}
    table = table.filter(feats & set(table.ids(axis='observation')), axis='observation')
    if depth is not None:
        table = table.filter(lambda v, i, m: v.sum() >= depth).remove_empty()
    return table



@click.command()
@click.option('--taxonomy', type=click.Path(exists=True), required=True)
@click.option('--tree', type=click.Path(exists=True), required=True)
@click.option('--table-16s', type=click.Path(exists=True), required=True)
@click.option('--table-wgs', type=click.Path(exists=True), required=True)
@click.option('--metadata', type=click.Path(exists=True), required=True)
@click.option('--output', type=click.Path(exists=False), required=True)
@click.option('--threads', type=int, required=True)
@click.option('--depth-16S', type=int, default=10000)
@click.option('--depth-WGS', type=int, default=1000000)
def process(taxonomy, tree, table_16s, table_wgs, output, threads, metadata,
            depth_16s, depth_wgs):
    taxonomy = qiime2.Artifact.load(taxonomy).view(pd.DataFrame)
    tree_ar = qiime2.Artifact.load(tree)
    tree = tree_ar.view(skbio.TreeNode)
    table_16s_ar = qiime2.Artifact.load(table_16s)
    table_wgs_ar = qiime2.Artifact.load(table_wgs)
    table_16s = table_16s_ar.view(biom.Table)
    table_wgs = table_wgs_ar.view(biom.Table)
    metadata = qiime2.Metadata.load(metadata).to_dataframe()

    #table_16s = subsample(table_16s, tree, depth_16s)
    #table_wgs = subsample(table_wgs, tree, depth_wgs)
    table_16s = filter_table(table_16s, tree, depth_16s)
    table_wgs = filter_table(table_wgs, tree, depth_wgs)
    table, metadata = concat(table_16s, table_wgs, metadata)

    metadata.to_csv(output + '.metadata.tsv', sep='\t', index=True, header=True)

    table_ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', table)
    table_ar.save(output + '.feature_table.qza')

    taxonomy['genus'] = taxonomy['Taxon'].apply(lambda x: x.split('; ')[-2])
    genus_table = as_genus(table, taxonomy)
    genus_table_ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', genus_table)

    
    wu_feature_dm, = diversity.actions.beta_phylogenetic(table_ar, tree_ar,
                                                         threads=threads,
                                                         metric='weighted_normalized_unifrac')
    wu_feature_pc, = diversity.actions.pcoa(wu_feature_dm, number_of_dimensions=5)
    wu_feature_emp, = emperor.actions.plot(wu_feature_pc, qiime2.Metadata(metadata))
    wu_feature_dm.save(output + '.asv.weighted.qza')
    wu_feature_pc.save(output + '.asv.weighted.pc.qza')
    wu_feature_emp.save(output + '.asv.weighted.pc.emp.qzv')


if __name__ == '__main__':
    process()
