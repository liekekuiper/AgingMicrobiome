#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 164G
#SBATCH -J bowtie2_build
#SBATCH -o Logs/slurm-%x-%j.out
#SBATCH -e Logs/slurm-%x-%j.err
#### One-time setup, adapted from PreProcessingShotgun/make_bowtie2_db.sbatch:
#### builds the bowtie2 index from the WoL2 reference genomes
#### ($WOL2DB/all.fna) into the location config.sh's DB variable points at.
####
#### Deliberately NOT wired into submit_batches.sh's automatic job chain
#### (unlike the original pipeline's single wrapper script, which rebuilt
#### the index on every run): the index only depends on the WoL2 reference,
#### not on your sample manifest, so rebuilding it for every new batch of
#### cohort samples would be wasted 164G/16-core compute. Run this once,
#### manually, before your first submit_batches.sh run -- and again only if
#### you update the WoL2 reference itself.
####
#### Log path is hardcoded to Logs/ above (matching config.sh's LOGS_ROOT
#### default) rather than read from config.sh, because SLURM needs the
#### #SBATCH -o/-e paths at submission time, before the script (and its
#### `source config.sh`) ever runs. That directory must already exist when
#### you submit -- unlike submit_batches.sh, which creates it for you, run:
####   mkdir -p Logs
#### first if it doesn't exist yet (or adjust the two lines above to match
#### a non-default LOGS_ROOT).
####
#### Prerequisite: download and decompress the reference first (see the
#### README):
####   mkdir -p wol2 && cd wol2
####   curl -L -O https://ftp.microbio.me/pub/wol2/genomes/all.fna.xz
####   unxz all.fna.xz
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
cd "$SLURM_SUBMIT_DIR"

WOL2DB_ABS=$(readlink -m "$WOL2DB")
DB_ABS=$(readlink -m "$DB")
mkdir -p "$(dirname "$DB_ABS")"

if [[ ! -f "$WOL2DB_ABS/all.fna" ]]; then
  echo "ERROR: $WOL2DB_ABS/all.fna not found." >&2
  echo "Download and decompress it first:" >&2
  echo "  cd $WOL2DB_ABS && curl -L -O https://ftp.microbio.me/pub/wol2/genomes/all.fna.xz && unxz all.fna.xz" >&2
  exit 1
fi

echo "# Building bowtie2 index from $WOL2DB_ABS/all.fna -> $DB_ABS"
echo "# "$(date)
bowtie2-build --threads 16 "$WOL2DB_ABS/all.fna" "$DB_ABS"
echo "# Done. DB= in config.sh should point to: $DB"
echo "# "$(date)
