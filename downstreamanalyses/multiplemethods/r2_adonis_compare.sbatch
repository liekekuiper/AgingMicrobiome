#!/bin/bash
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --time=8-00:00:00
#SBATCH -J pcoa
#SBATCH --mem=64G

source activate qiime2-2023.7

qiime diversity adonis \
--i-distance-matrix MrOS.asv.weighted.qza \
--m-metadata-file MrOS.metadata_with_sample.tsv \
--o-visualization weighted_prep_significance.qzv \
--p-formula preparation

qiime diversity adonis \
--i-distance-matrix MrOS.asv.weighted.qza \
--m-metadata-file MrOS.metadata_with_sample.tsv \
--o-visualization weighted_sample_significance.qzv \
--p-formula sample

