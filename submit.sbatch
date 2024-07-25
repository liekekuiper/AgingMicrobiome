#!/bin/bash
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --time=8:00:00
#SBATCH -J mb_aging_downstream
#SBATCH --mem=64G

# Activating environment
source activate qiime2-2023.7


python fixed_analysiscode2.py \
         --tree=${tree} \
         --taxonomy=${taxonomy} \
         --metadata=${metadata} \
         --feature_table=${feature_table} \
         --threads=${threads} \
         --modsfile=${modsfile} \
         --outsfile=${outsfile} \
         --cohortname=${cohortname} \
         --subsfile=${subsfile}
