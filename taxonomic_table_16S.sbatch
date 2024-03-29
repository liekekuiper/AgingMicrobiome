#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time 24:00:00
#SBATCH --mem 64gb
#SBATCH -J  taxonomy_from_table

source activate qiime2-2023.7

set -x
set -e

# Get taxonomy based on backbone

#qiime 
qiime greengenes2 taxonomy-from-table \
     --i-reference-taxonomy 2022.10.taxonomy.asv.nwk.qza \
     --i-table feature.table.gg2-2022.10.qza \
     --o-classification df.gg2.taxonomy.qza
