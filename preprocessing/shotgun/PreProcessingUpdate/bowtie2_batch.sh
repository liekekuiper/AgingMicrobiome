#!/bin/bash
#SBATCH -c 16
#SBATCH --mem 128G
#SBATCH -J bowtie2_batch
#SBATCH -o slurm-%x-%A_%a.out
#SBATCH -e slurm-%x-%A_%a.err

#### This file uses the manifest file to create batches to run woltka.sbatch on
#### woltka.sbatch needs .sam files which are too large for most HPC-clusters
#### Hence, we split them in batches and make them into .bam files after they are run
# What it does:
# Reads the manifest file and runs bowtie2 for each sample in the batch.
# Produces .sam files in the batch align directory.
# .sam files are converted to .bam by sam2bam_batch.sh after woltka completes.

set -euo pipefail

MANIFEST=manifest
BATCH_SIZE=10
DB=/path/to/bowtie2/db         # Change to your bowtie2 database index path
BATCH_ID=$SLURM_ARRAY_TASK_ID  # 1-based batch index passed via job array
BATCH_DIR=$(printf "Batch%02d" "$BATCH_ID")
ALIGN="$BATCH_DIR/align"

cd "$SLURM_SUBMIT_DIR"
mkdir -p "$ALIGN"

START=$(( (BATCH_ID - 1) * BATCH_SIZE + 1 ))
END=$(( START + BATCH_SIZE - 1 ))

mapfile -t ROWS < <(tail -n +2 "$MANIFEST")

for ((i=START; i<=END && i<=${#ROWS[@]}; i++)); do
  ROW="${ROWS[$((i - 1))]}"
  SAMPLE=$(cut -f1 <<< "$ROW")
  READS=$(cut -f2 <<< "$ROW")
  OUT="$ALIGN/$SAMPLE.sam"
  [ -f "$OUT" ] && continue    # Skip if output already exists
  bowtie2 -x "$DB" -U "$READS" -S "$OUT" -p 16
done
