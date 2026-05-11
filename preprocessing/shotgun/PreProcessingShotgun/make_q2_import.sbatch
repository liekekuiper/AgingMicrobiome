#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time 8:00:00
#SBATCH --mem 64gb
#SBATCH -J q2_import

source activate qiime2-2023.7

set -x
set -e

qiime tools import \
    --input-path Output/output.biom \
    --output-path woltka.biom.qza \
    --type FeatureTable[Frequency]
