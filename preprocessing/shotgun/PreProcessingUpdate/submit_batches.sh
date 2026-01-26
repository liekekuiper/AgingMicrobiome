#!/bin/bash
set -euo pipefail

MANIFEST=manifest
BATCH_SIZE=10

TOTAL=$(($(wc -l < "$MANIFEST") - 1))
NBATCH=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))

for ((b=1; b<=NBATCH; b++)); do
  echo "Submitting batch $b"

  BOWTIE_JOB=$(sbatch --array=$b bowtie2_batch.sh | awk '{print $4}')
  WOLTKA_JOB=$(sbatch --dependency=afterok:$BOWTIE_JOB woltka_batch.sh $b | awk '{print $4}')
  sbatch --dependency=afterok:$WOLTKA_JOB sam2bam_batch.sh $b
done
