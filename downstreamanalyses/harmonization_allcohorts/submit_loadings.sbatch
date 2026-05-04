#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time 8:00:00
#SBATCH -J beta
#SBATCH --mem 64g

source activate qiime2-2023.7

python create_loadings_top_genera.py \
    --taxonomy ./MrOS/MrOS_Shotgun_Sequencing/2022.10.taxonomy.asv.tsv.qza \
    --tables ./FHS/16S/feature.table.gg2-2022.10.qza \
    --tables ./FHS/metagenomics/reindexed-table.qza \
    --tables ./MrOS/16S/renamed_feature.table.gg2-2022.10.qza \
    --tables ./MrOS/MrOS_Shotgun_Sequencing/renamed_feature.table.gg2-2022.10.qza \
    --tables ./Analyses/feature.table.gg2-2022.10.qza \
     --tables ./SOL/metagenomics/filtered_renamed_shot_feature.table.gg2-2022.10.qza \
    --tables ./SOL/16S/filtered_renamed_16s_feature.table.gg2-2022.10.qza \
    --cohort-names FHS16S \
    --cohort-names FHSshotgun \
    --cohort-names MrOS16S \
    --cohort-names MrOSshotgun \
    --cohort-names RS \
    --cohort-names SOLshotgun \
    --cohort-names SOL16S \
    --output-dir ./ \
    --top-n 20
