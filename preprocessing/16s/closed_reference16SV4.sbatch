#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 64gb
#SBATCH -J closedref

source activate qiime2-2023.7

set -x
set -e

qiime greengenes2 filter-features \
    --i-feature-table data.deblur.150nt.qza \
    --i-reference ../2022.10.taxonomy.asv.nwk.qza \
    --o-filtered-feature-table feature.table.gg2-2022.10.qza 
