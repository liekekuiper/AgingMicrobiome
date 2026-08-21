#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time 8:00:00
#SBATCH --mem 64G
#SBATCH -J q2_import
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err
#### (submit_batches.sh overrides these two via --output/--error so logs
#### land in $LOGS_ROOT; only used as-is if you sbatch this file directly.)
#### Adapted from PreProcessingShotgun/make_q2_import.sbatch: imports the
#### final merged OGU table into QIIME2 as a FeatureTable[Frequency]
#### artifact. Where the original imported the single-batch Output/output.biom,
#### this imports Merged/Output/output.biom -- the table merge_batches.sh
#### produces across the whole cohort.
####
#### Runs after merge_batches.sh finishes (submit_batches.sh wires up this
#### dependency automatically).
####
#### The GreenGenes2 follow-up step (../closed_referenceShotgun.sbatch)
#### expects woltka.biom.qza in the main pipeline folder (one level above
#### Scripts/, alongside Data/, wol2/, db/), not buried under
#### woltka_files/Merged/Output/ -- so once the .qza is built, this script
#### also copies it up there.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
cd "$SLURM_SUBMIT_DIR"
source activate "$QIIME2_ENV"

OUTPUT_ROOT_ABS=$(readlink -m "$OUTPUT_ROOT")
IN="$OUTPUT_ROOT_ABS/Merged/Output/output.biom"
OUT="$OUTPUT_ROOT_ABS/Merged/Output/woltka.biom.qza"
MAIN_FOLDER=$(readlink -m "$SLURM_SUBMIT_DIR/..")

if [[ ! -f "$IN" ]]; then
  echo "ERROR: $IN not found -- did merge_batches.sh run and finish successfully?" >&2
  exit 1
fi

echo "# Importing $IN -> $OUT"
echo "# "$(date)
qiime tools import \
  --input-path "$IN" \
  --output-path "$OUT" \
  --type FeatureTable[Frequency]

cp "$OUT" "$MAIN_FOLDER/woltka.biom.qza"
echo "# Also copied to $MAIN_FOLDER/woltka.biom.qza -- run closed_referenceShotgun.sbatch from there next."
echo "# Done."
echo "# "$(date)
