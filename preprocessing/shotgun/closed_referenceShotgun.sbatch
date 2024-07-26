#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time 8:00:00
#SBATCH --mem 64G
#SBATCH -J closed_ref

source activate qiime2-2023.7

set -x
set -e

qiime greengenes2 filter-features \
    --i-feature-table woltka.biom.qza \
    --i-reference 2022.10.taxonomy.asv.nwk.qza \
    --o-filtered-feature-table feature.table.gg2-2022.10.qza 
