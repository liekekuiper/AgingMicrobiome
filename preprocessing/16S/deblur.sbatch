#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 64gb
#SBATCH -J deblur

source activate qiime2-2023.7

set -x
set -e

qiime deblur denoise-16S \
    --i-demultiplexed-seqs data.seqs.read1.qza \
    --p-no-hashed-feature-ids \
    --p-jobs-to-start 16 \
    --p-min-reads 1 \
    --p-trim-length 150 \
    --o-table data.deblur.150nt.qza \
    --o-representative-sequences data.repset.150nt.qza \
    --o-stats data.deblur.stats.150nt.qza
