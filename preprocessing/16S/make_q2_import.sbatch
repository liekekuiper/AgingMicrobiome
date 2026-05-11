#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time 8:00:00
#SBATCH --mem 8gb
#SBATCH -J import

source activate qiime2-2023.7

set -x
set -e

qiime tools import \
    --type SampleData[SequencesWithQuality] \
    --input-path manifest \
    --output-path data.seqs.read1.qza \
    --input-format SingleEndFastqManifestPhred33V2
