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
metadata=t #path to metadata
feature_table= #path to feature table
tax_table= #path to taxonomy table
threads= #number of threads
modsfile= #path to mods_agingmicrobiome.txt
outsfile= #path to  outs_agingmicrobiome.txt
cohortname= #path to cohort.txt
subsfile= #path to subsets_agingmicrobiome.txt

# Submitting job
sbatch --export tree=${tree},taxonomy=${taxonomy},metadata=${metadata},feature_table=${feature_table},tax_table=${tax_table},threads=${threads},modsfile=${modsfile},outsfile=${outsfile},cohortname=${cohortname},subsfile=${subsfile},label=${label} submit.sbatch
