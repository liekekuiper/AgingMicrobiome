#!/bin/bash
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -J woltka_batch
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err


### If you want functions run the following two lines (without the #) in your wol2 folder
# wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/function/
# wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/genomes/

#if not turn on --no-fun

set -euo pipefail

BATCH_ID=$1
BATCH_DIR=$(printf "Batch%02d" "$BATCH_ID")

cd "$SLURM_SUBMIT_DIR/$BATCH_DIR"

./wol2sop.sh \
  -d ../wol2 \ #Change into your wol2 directory
  -i align \
  -o Output \
  -f biom \
  --no-tax
#--no-fun
