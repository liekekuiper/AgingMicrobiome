#!/bin/bash
#SBATCH -c 4
#SBATCH --mem 126G
#SBATCH -J merge_batches
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err
#### (submit_batches.sh overrides these two via --output/--error so logs
#### land in $LOGS_ROOT; only used as-is if you sbatch this file directly.)
#### Memory note: merging orf.biom across all batches (ORF-level, no
#### collapsing) produces a huge matrix -- at 1347 samples this was already
#### 3.5M+ features and OOM-killed at 32G. 126G is a starting guess; check
#### your cluster's node memory limits and adjust, and watch `seff <jobid>`
#### / sacct after a run to see actual peak (MaxRSS) and tune from there.
#### The other tables (uniref/go/pfam/kegg/metacyc/eggnog, all collapsed to
#### much smaller feature spaces) are far cheaper -- orf.biom is the one
#### driving this requirement.
#### align_classify_batch.sh already combines every sample WITHIN a batch
#### into one shared table (all samples in that batch become columns of the
#### same orf.biom, uniref.biom, etc. -- nothing is per-sample). But each
#### batch still writes to its own BatchNN/Output/ directory, so across
#### batches you end up with several separate tables, each covering only
#### that batch's samples.
####
#### This script closes that gap: it walks the UNION of relative file paths
#### found across every BatchNN/Output/ (orf.biom, uniref/uniref.biom,
#### go/all.biom, kegg/pathway.biom, ... all ~68 tables produced by the
#### functional cascade -- not just whatever happens to be in batch 1,
#### since a batch that failed partway through would otherwise silently
#### shrink what gets merged) and merges each path across every batch that
#### has it into one final table under Merged/Output/, using `woltka
#### merge`. Because each batch has a disjoint set of samples, the merge is
#### a union: no sample data is summed or altered, you simply end up with
#### one table per profile type that has every sample in the manifest as a
#### column.
####
#### Run this after all align_classify_batch.sh array tasks have finished
#### (submit_batches.sh wires up this dependency automatically).
####
#### Resumable: a table already present under Merged/Output/ is skipped on
#### rerun, not recomputed. This is what makes it cheap to flip MERGE_ORF
#### on later and re-sbatch this script directly (without going through
#### submit_batches.sh/run_pipeline.sh, and without re-running any
#### batches) -- only the newly-enabled orf.biom gets merged, everything
#### already sitting in Merged/Output/ is left untouched. Set
#### FORCE_MERGE=1 in config.sh if you ever need to force a full re-merge
#### of everything (e.g. after fixing and re-running one specific batch).
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
cd "$SLURM_SUBMIT_DIR"

# Same OUTPUT_ROOT as align_classify_batch.sh, so this looks in the same
# place those jobs wrote their BatchNN/Output/ folders to.
OUTPUT_ROOT_ABS=$(readlink -m "$OUTPUT_ROOT")

TOTAL=$(($(wc -l < "$MANIFEST") - 1))
NBATCH=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))

# Build the UNION of relative output paths across ALL batches (not just
# batch 1 -- a batch that failed partway through the cascade would
# otherwise silently limit what gets merged).
declare -A SEEN
FOUND_ANY_BATCH=
for ((b=1; b<=NBATCH; b++)); do
  BD="$OUTPUT_ROOT_ABS/$(printf "Batch%02d" "$b")/Output"
  [[ -d "$BD" ]] || continue
  FOUND_ANY_BATCH=1
  while IFS= read -r -d '' f; do
    SEEN["${f#"$BD"/}"]=1
  done < <(find "$BD" -type f -print0)
done

if [[ -z "$FOUND_ANY_BATCH" ]]; then
  echo "ERROR: no BatchNN/Output directories found under $OUTPUT_ROOT_ABS -- did the batches run and finish successfully?" >&2
  exit 1
fi

# orf.biom is the raw, uncollapsed per-ORF table -- by far the heaviest
# thing to merge (it's what OOM-killed this job at 32G). Skipped by default
# (MERGE_ORF unset in config.sh); the per-batch BatchNN/Output/orf.biom
# files are untouched either way, only the merged copy is skipped.
if [[ -z "${MERGE_ORF:-}" ]] && [[ -n "${SEEN[orf.biom]:-}" ]]; then
  unset 'SEEN[orf.biom]'
  echo "# Skipping merge of orf.biom (MERGE_ORF is off in config.sh -- per-batch orf.biom files are untouched)."
fi

MERGED="$OUTPUT_ROOT_ABS/Merged/Output"
mkdir -p "$MERGED"

echo "# Merging $NBATCH batches into $MERGED (${#SEEN[@]} distinct output tables found)"
echo "# "$(date)

# NOTE: each table is merged independently below, and a failure on one
# does NOT abort the rest (`if ! cmd; then ...; fi` doesn't trigger
# `set -e`) -- otherwise a single bad merge (e.g. orf.biom) would silently
# kill every table queued after it, leaving Merged/Output only partially
# populated with no explanation.
FAILED=()
SKIPPED=()
mapfile -t RELS < <(printf '%s\n' "${!SEEN[@]}" | sort)
for rel in "${RELS[@]}"; do
  if [[ -z "${FORCE_MERGE:-}" ]] && [[ -f "$MERGED/$rel" ]]; then
    SKIPPED+=("$rel")
    continue
  fi

  paths=()
  for ((b=1; b<=NBATCH; b++)); do
    p="$OUTPUT_ROOT_ABS/$(printf "Batch%02d" "$b")/Output/$rel"
    [[ -f "$p" ]] && paths+=("$p")
  done
  [[ ${#paths[@]} -eq 0 ]] && continue

  mkdir -p "$(dirname "$MERGED/$rel")"

  if [[ ${#paths[@]} -eq 1 ]]; then
    # Only one batch produced this file (e.g. a single-batch run) -- nothing
    # to merge, just place it.
    if ! cp "${paths[0]}" "$MERGED/$rel"; then
      echo "WARNING: failed to copy $rel -- skipping." >&2
      FAILED+=("$rel")
      continue
    fi
  else
    args=()
    for p in "${paths[@]}"; do args+=(-i "$p"); done
    echo "# Merging $rel across ${#paths[@]} batches..."
    if ! woltka merge "${args[@]}" -o "$MERGED/$rel"; then
      echo "WARNING: woltka merge failed for $rel -- skipping (see error above)." >&2
      FAILED+=("$rel")
      continue
    fi
  fi

  if [[ ${#paths[@]} -lt $NBATCH ]]; then
    echo "# NOTE: $rel was only found in ${#paths[@]}/$NBATCH batches -- the missing batch(es) may have failed partway through. Check their slurm .err logs."
  fi
done

DONE=$(( ${#RELS[@]} - ${#FAILED[@]} - ${#SKIPPED[@]} ))
echo "# Merge completed: $DONE table(s) merged, ${#SKIPPED[@]} already up to date (skipped), out of ${#RELS[@]} total, into $MERGED/"
if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "# ${#FAILED[@]} table(s) FAILED to merge -- see WARNING lines above for the underlying error:" >&2
  printf '#   %s\n' "${FAILED[@]}" >&2
fi
echo "# "$(date)
