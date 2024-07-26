#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 64gb
#SBATCH -J closedref

source activate qiime2-2023.7

set -x
set -e

python get_rep_set.py data.deblur.150nt.qza data.deblur.150nt.repset.qza
qiime greengenes2 non-v4-16s \
    --i-table data.deblur.150nt.qza \
    --i-sequences data.deblur.150nt.repset.qza \
    --i-backbone ../2022.10.backbone.full-length.fna.qza \
    --p-threads 16 \
    --o-mapped-table feature.table.gg2-2022.10.qza \
    --o-representatives df.gg2-2022.10-cref99.repset.qza
