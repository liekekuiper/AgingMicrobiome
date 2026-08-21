#!/bin/bash
#### Submits the align+classify batches, then merge_batches.sh, then
#### (if RUN_Q2_IMPORT is on in config.sh) make_q2_import.sh -- everything
#### EXCEPT the one-time bowtie2 index build (see make_bowtie2_db.sh). Run
#### this directly for reruns once your index already exists; for a full
#### first-time run (index included), use `bash run_pipeline.sh` instead,
#### which calls this script after submitting make_bowtie2_db.sh.
####
#### EXTRA_DEPENDENCY (optional env var): if set, e.g.
#### "afterok:12345", this is added as a dependency on the very first
#### batch of align_classify_batch.sh jobs -- this is how run_pipeline.sh
#### makes the batches wait for the index build without this script
#### needing to know anything about that step itself.
set -euo pipefail
source config.sh
TOTAL=$(($(wc -l < "$MANIFEST") - 1))
NBATCH=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))

# Every script below has its own #SBATCH -o/-e as a fallback (in case you
# ever sbatch one directly, e.g. to retry a single failed batch), but here
# we override those via CLI flags so all logs from a submit_batches.sh run
# land in one place: $LOGS_ROOT. --output/--error on the sbatch command
# line take precedence over the #SBATCH header lines in the script.
LOGS_ABS=$(readlink -m "$LOGS_ROOT")
mkdir -p "$LOGS_ABS"

BATCH_DEPS=()
[[ -n "${EXTRA_DEPENDENCY:-}" ]] && BATCH_DEPS+=("--dependency=$EXTRA_DEPENDENCY")

JOB_IDS=()
for ((b=1; b<=NBATCH; b++)); do
  echo "Submitting batch $b"
  JOB_ID=$(sbatch --parsable --array=$b "${BATCH_DEPS[@]}" \
    --output="$LOGS_ABS/slurm-%x-%A_%a.out" \
    --error="$LOGS_ABS/slurm-%x-%A_%a.err" \
    align_classify_batch.sh)
  JOB_IDS+=("$JOB_ID")
done

# One job merges every batch's output tables into a single set of tables
# covering all samples, once every batch has finished successfully.
DEPENDENCY="afterok:$(IFS=:; echo "${JOB_IDS[*]}")"
echo "Submitting merge_batches.sh (depends on: $DEPENDENCY)"
MERGE_JOB_ID=$(sbatch --parsable --dependency="$DEPENDENCY" \
  --output="$LOGS_ABS/slurm-%x-%j.out" \
  --error="$LOGS_ABS/slurm-%x-%j.err" \
  merge_batches.sh)

# Finally, import the merged OGU table into QIIME2, once the merge is
# done -- unless RUN_Q2_IMPORT is off in config.sh, in which case we stop
# here and no .qza is produced (or needed: QIIME2_ENV is never touched).
if [[ -n "${RUN_Q2_IMPORT:-}" ]]; then
  echo "Submitting make_q2_import.sh (depends on: afterok:$MERGE_JOB_ID)"
  sbatch --parsable --dependency="afterok:$MERGE_JOB_ID" \
    --output="$LOGS_ABS/slurm-%x-%j.out" \
    --error="$LOGS_ABS/slurm-%x-%j.err" \
    make_q2_import.sh
else
  echo "RUN_Q2_IMPORT is off in config.sh -- skipping make_q2_import.sh."
fi
