#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -J  bowtie2

# Create the align directory if it doesn't exist
mkdir -p align

# Iterate through each file in the Data directory
for file in Data/*; do
    if [ -f "$file" ]; then
        # Get the file name without the extension
        filename=$(basename -- "$file")
        filename_no_ext="${filename%.*}"

        # Run bowtie2 on the file and save the output .sam file in the align directory
        bowtie2 -x db  -U "$file" -S "align/${filename_no_ext}.sam"
    fi
done

