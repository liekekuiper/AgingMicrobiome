#!/bin/bash

#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 164G
#SBATCH -J bowtie2_build

bowtie2-build --threads 9 ./wol2/all.fna db
