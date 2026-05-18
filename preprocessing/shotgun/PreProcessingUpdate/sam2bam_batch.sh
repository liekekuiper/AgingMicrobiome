#!/bin/bash
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -J sam2bam_batch
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

# Converts .sam files to sorted, indexed .bam files and removes the originals.
# Runs after woltka_batch.sh completes to free up disk space.

set -euo pipefail

BATCH_ID=$1
BATCH_DIR=$(printf "Batch%02d" "$BATCH_ID")
ALIGN="$SLURM_SUBMIT_DIR/$BATCH_DIR/align"

shopt -s nullglob
sams=("$ALIGN"/*.sam)

if [[ ${#sams[@]} -eq 0 ]]; then
  echo "No .sam files found in $ALIGN, nothing to do."
  exit 0
fi

for sam in "${sams[@]}"; do
  bam="${sam%.sam}.bam"
  samtools view -@ 16 -b "$sam" \
    | samtools sort -@ 16 -o "$bam"
  samtools index "$bam"
  rm -f "$sam"    # Remove .sam after successful conversion to save disk space
done
