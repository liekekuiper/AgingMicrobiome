#!/bin/bash

# Setting up variables
version=2022.10
taxonomy=${version}.taxonomy.asv.tsv.qza
tree=${version}.phylogeny.asv.nwk.qza

# Check and download the taxonomy file if it doesn't exist
if [ ! -f "$taxonomy" ]; then
  echo "Downloading $taxonomy..."
  curl -O "http://ftp.microbio.me/greengenes_release/${version}/${taxonomy}"
else
  echo "$taxonomy already exists. Skipping download."
fi

# Check and download the tree file if it doesn't exist
if [ ! -f "$tree" ]; then
  echo "Downloading $tree..."
  curl -O "http://ftp.microbio.me/greengenes_release/${version}/${tree}"
else
  echo "$tree already exists. Skipping download."
fi

label= #Fill in 16s or shotgun
metadata=../metadata.txt #path to metadata
feature_table=feature.table.gg2-2022.10.qza #path to feature table
tax_table=df.gg2.taxonomy.qza #path to taxonomy table
threads=6 #number of threads
modsfile=../mods_agingmicrobiome.txt #path to mods_agingmicrobiome.txt
outsfile=../outs_agingmicrobiome.txt #path to  outs_agingmicrobiome.txt
cohortname=../cohort.txt #path to cohort.txt
subsfile=../subsets_agingmicrobiome.txt #path to subsets_agingmicrobiome.txt
factors="NA" #list study-specific categorical variables as "factor1,factor2" if none specify "NA" 

# Submitting job
sbatch --export tree=${tree},taxonomy=${taxonomy},metadata=${metadata},feature_table=${feature_table},tax_table=${tax_table},threads=${threads},modsfile=${modsfile},outsfile=${outsfile},cohortname=${cohortname},subsfile=${subsfile},label=${label},factors="'${factors}'" submit.sbatch
