#!/bin/bash
#### Runs the full PreProcessingShotgunUpdate pipeline end-to-end with one
#### command -- adapted from PreProcessingShotgun's original wrapper:
####   j=$(sbatch --parsable make_bowtie2_db.sbatch)
####   j=$(sbatch --parsable --dependency afterok:${j} bowtie2.sbatch)
####   j=$(sbatch --parsable --dependency afterok:${j} woltka.sbatch)
####   j=$(sbatch --parsable --dependency afterok:${j} make_q2_import.sbatch)
####
#### Same idea, adapted to this pipeline's streaming batch/merge design:
####   1. make_bowtie2_db.sh          (skipped if the index already exists,
####                                   see below -- the original always
####                                   rebuilt it, which is wasteful once
####                                   you already have one)
####   2. align_classify_batch.sh, once per batch  \
####   3. merge_batches.sh                          } via submit_batches.sh
####   4. make_q2_import.sh                        /
####
#### Run from inside Scripts/:
####   bash run_pipeline.sh
####
#### To force a fresh index rebuild even if one already exists, remove (or
#### rename) the existing $DB.*.bt2/.bt2l files first, or just run
#### `sbatch make_bowtie2_db.sh` yourself before this script.
set -euo pipefail
source config.sh
cd "$(dirname "${BASH_SOURCE[0]}")"

LOGS_ABS=$(readlink -m "$LOGS_ROOT")
mkdir -p "$LOGS_ABS"

DB_ABS=$(readlink -m "$DB")

# --- 1. bowtie2 index (one-time; skipped if it already exists) ---
if [[ -f "$DB_ABS.1.bt2" || -f "$DB_ABS.1.bt2l" ]]; then
  echo "# Bowtie2 index already exists at ${DB_ABS}.1.bt2[l] -- skipping make_bowtie2_db.sh."
  EXTRA_DEPENDENCY=""
else
  echo "# Submitting make_bowtie2_db.sh"
  DB_JOB_ID=$(sbatch --parsable \
    --output="$LOGS_ABS/slurm-%x-%j.out" \
    --error="$LOGS_ABS/slurm-%x-%j.err" \
    make_bowtie2_db.sh)
  EXTRA_DEPENDENCY="afterok:$DB_JOB_ID"
  echo "# make_bowtie2_db.sh submitted as job $DB_JOB_ID; remaining steps will wait for it."
fi

# --- 2-4. align+classify per batch, merge, QIIME2 import ---
export EXTRA_DEPENDENCY
bash submit_batches.sh
