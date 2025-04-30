#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 64gb
#SBATCH -J abundance

source activate qiime2-2023.7

set -x
set -e

qiime taxa collapse \
    --i-table feature.table.gg2-2022.10.qza \
    --i-taxonomy df.gg2.taxonomy.qza \
    --p-level 7 \
    --o-collapsed-table species-table.qza 

qiime feature-table relative-frequency \
    --i-table species-table.qza \
    --o-relative-frequency-table rel-species-table.qza

python rel_ab_prev.py
