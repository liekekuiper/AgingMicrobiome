#!/bin/bash
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -J woltka_batch
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

# Runs wol2sop.sh on the aligned .sam files for a given batch.
# To include functional profiling, download the following into your wol2 directory:
#   wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/function/
#   wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/genomes/
#   wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/proteins/
# If not interested in functions remove # in front of --no-fun below.

set -euo pipefail

BATCH_ID=$1
BATCH_DIR=$(printf "Batch%02d" "$BATCH_ID")

"$SLURM_SUBMIT_DIR/wol2sop.sh" \
  -d "$SLURM_SUBMIT_DIR/wol2" \
  -i "$SLURM_SUBMIT_DIR/$BATCH_DIR/align" \
  -o "$SLURM_SUBMIT_DIR/$BATCH_DIR/Output" \
  -f biom \
  --no-tax
# --no-fun
