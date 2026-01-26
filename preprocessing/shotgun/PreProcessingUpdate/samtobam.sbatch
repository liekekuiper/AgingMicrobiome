#!/bin/bash
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -J sam2bam_batch
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

set -euo pipefail

BATCH_ID=$1
BATCH_DIR=$(printf "Batch%02d" "$BATCH_ID")
ALIGN="$SLURM_SUBMIT_DIR/$BATCH_DIR/align"

for sam in "$ALIGN"/*.sam; do
  [ -e "$sam" ] || exit 0
  bam="${sam%.sam}.bam"

  samtools view -@ 16 -bS "$sam" \
    | samtools sort -@ 16 -o "$bam"

  samtools index "$bam"
  rm -f "$sam"
done
