# PreProcessingShotgunUpdate

Faster, streaming version of the shotgun metagenomics preprocessing
pipeline (Bowtie2 + Woltka). See `../README.md` (the `shotgun/` directory
README) for how this relates to the original `PreProcessingShotgun/`
pipeline used in the manuscript, and which one to use.

This pipeline expects already quality-controlled (kneaddata-processed)
reads as input -- QC is assumed to be done separately, before this
pipeline.

## 1. Install dependencies

```
pip install woltka
```

You'll also need `bowtie2` and `gzip`/`pigz` available (as modules or on
your `$PATH`).

## 2. Set up the WoL2 database

Download the Web of Life reference genomes, needed to build your bowtie2
index:

```
mkdir -p wol2
cd wol2
curl -L -O https://ftp.microbio.me/pub/wol2/genomes/all.fna.xz
unxz all.fna.xz
```

Set `DB=` in `config.sh` to where you want the bowtie2 index to live. You
don't need to build it yourself: `run_pipeline.sh` (its own step 1) builds
it automatically the first time, via `Scripts/make_bowtie2_db.sh` (adapted
from `PreProcessingShotgun/make_bowtie2_db.sbatch`), and skips that step
on later runs once the index already exists. Only run it manually if you
want to (re)build the index on its own -- e.g. after updating the WoL2
reference:

```
cd Scripts
mkdir -p Logs   # make_bowtie2_db.sh can't create this itself in time, see its header
sbatch make_bowtie2_db.sh
```

Then download the function and protein reference data used for
classification:

```
wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/function/
wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/proteins/
```

Note: you do **not** need the rest of `genomes/` beyond `all.fna` above --
this pipeline doesn't use taxonomic-rank or genome-size normalization, so
`genomes/length.map` and friends aren't required.

## 3. Build the sample manifest

Make a file `manifest` with a column `sample-id` and a column
`absolute-filepath`, containing the sample IDs and the paths to each
sample's QC'd (kneaddata) forward reads: a header line, then one
tab-separated `sample-id`/`absolute-filepath` pair per line.

Sample IDs must not contain `_`, since this pipeline demultiplexes samples
by splitting on the first underscore in each read ID -- an underscore in
a sample ID would silently misattribute reads to the wrong sample.

## 4. Configure and run

Suggested layout:

```
PreProcessingShotgunUpdate/
├── Scripts/
│   ├── config.sh
│   ├── run_pipeline.sh
│   ├── make_bowtie2_db.sh
│   ├── submit_batches.sh
│   ├── align_classify_batch.sh
│   ├── merge_batches.sh
│   └── make_q2_import.sh
├── Data/            # QC'd .fq.gz read files
├── db/              # bowtie2 index, built in step 2
├── wol2/            # WoL2 database, downloaded in step 2
├── woltka_files/    # batch + merged output (OUTPUT_ROOT, step 4) -- many, large files
└── Logs/            # SLURM stdout/stderr logs (LOGS_ROOT, step 4)
```

Edit `Scripts/config.sh`:

- `MANIFEST`, `BATCH_SIZE`, `DB`, `WOL2DB` -- paths and batch size for your
  run.
- `OUTPUT_ROOT` -- where batch and merged output get written. Defaults to
  `woltka_files` (recommended -- this pipeline produces many, large output
  files per batch, best kept out of `Scripts/` itself).
- `LOGS_ROOT` -- where SLURM stdout/stderr logs get written. Defaults to
  `Logs`; `submit_batches.sh` creates it automatically. If you run
  `make_bowtie2_db.sh` manually (step 2), create it yourself first
  (`mkdir -p Logs`) since that script's log path can't be read from
  `config.sh` until after the job has already started.
- `MERGE_ORF` -- leave off (default) unless you specifically need a single
  merged, per-ORF-resolution table across all samples; it's by far the
  most memory-hungry table to merge (a full cohort run needed ~300 GB just
  for that step).
- `RUN_Q2_IMPORT` -- on by default: imports the merged `output.biom` into
  QIIME2 as a `.qza` (`make_q2_import.sh`) after the merge finishes. Set
  to empty if you don't want a `.qza` -- that step is then skipped
  entirely and `QIIME2_ENV` below is never touched.
- `QIIME2_ENV` -- conda/mamba environment with QIIME2 installed, used by
  `make_q2_import.sh`. Unused if `RUN_Q2_IMPORT` is off.

Then run everything, from inside `Scripts/`, with one command:

```
cd Scripts
bash run_pipeline.sh
```

This mirrors `PreProcessingShotgun`'s original wrapper script (which
chained `make_bowtie2_db.sbatch` -> `bowtie2.sbatch` -> `woltka.sbatch` ->
`make_q2_import.sbatch` via `--dependency afterok`), adapted to this
pipeline's steps. `run_pipeline.sh` itself just handles the one-time index
build (skipping it if `$DB` already exists) and then calls
`submit_batches.sh`, which does the rest:

```
run_pipeline.sh
└─ make_bowtie2_db.sh          (skipped if $DB already exists)
   └─ submit_batches.sh
      ├─ align_classify_batch.sh, once per batch
      ├─ merge_batches.sh      (after all batches finish)
      └─ make_q2_import.sh     (after the merge finishes; skipped if
                                 RUN_Q2_IMPORT is off in config.sh)
```

If your index is already built and you just want to (re)run the batches
-- e.g. after adding new samples to the manifest -- you can call
`bash submit_batches.sh` directly instead and skip the index check.

## What this does

For each batch of samples (default 10), `align_classify_batch.sh`:

1. Aligns each sample with Bowtie2 and streams the combined output
   straight into Woltka -- no `.sam` file is ever written to disk.
2. Classifies OGU (genome-level) and ORF assignment concurrently, from
   that same stream.
3. Runs the full functional cascade on the resulting ORF table: UniRef,
   GO, Pfam, KEGG, MetaCyc, eggNOG, including external ID
   cross-references, KEGG/MetaCyc compound inference, and pathway/module
   coverage.

Once every batch finishes, `merge_batches.sh` merges every batch's tables
into one final set covering the whole cohort, written to `Merged/Output/`
under `OUTPUT_ROOT`.

Finally, if you want a `.qza` (`RUN_Q2_IMPORT=1` in `config.sh`, the
default), `make_q2_import.sh` imports `Merged/Output/output.biom` (the
merged OGU table) into QIIME2 as `Merged/Output/woltka.biom.qza`, a
`FeatureTable[Frequency]` artifact -- adapted from
`PreProcessingShotgun/make_q2_import.sbatch`, which did the same for the
single, non-batched `output.biom` the original pipeline produced. Set
`RUN_Q2_IMPORT=` (empty) if you don't need a `.qza` -- `Merged/Output/output.biom`
and the rest of the merged tables are produced either way; this only
controls the extra QIIME2 import step, and no QIIME2 environment is
needed if it's off.

## 5. Next step: GreenGenes2 harmonization

`Merged/Output/woltka.biom.qza` is **not** the final harmonized output --
GreenGenes2 calling is a separate follow-up step that lives one level up,
in the shared `shotgun/` directory (not duplicated inside
`PreProcessingShotgunUpdate/`), since it's the same step regardless of
whether you used this pipeline or `PreProcessingShotgun/`.

`make_q2_import.sh` already copies `woltka.biom.qza` up to the main
pipeline folder (alongside `Data/`, `wol2/`, `db/`) for you, so you don't
need to move anything yourself. From there, run **in this order**:

1. `../closed_referenceShotgun.sbatch` -- calls `woltka.biom.qza` against
   the GreenGenes2 backbone (`qiime greengenes2 filter-features`),
   producing `feature.table.gg2-2022.10.qza`.
2. `../taxonomic_table_Shotgun.sbatch` -- derives the corresponding
   taxonomy (`qiime greengenes2 taxonomy-from-table`) from that filtered
   table, producing `df.gg2.taxonomy.qza`.

These call GreenGenes2 release **2022.10**; see `../README.md` (the
`shotgun/` directory README) for how to switch to a newer release, and
for full GreenGenes2 details.

## What's deliberately different from `PreProcessingShotgun/`

`PreProcessingShotgun/` was run with `wol2sop.sh ... --no-tax --no-fun`:
OGU-level (genome) classification only, no ranked taxonomy, no ORF or
functional profiling at all, and as one single non-batched job for the
whole cohort.

- **ORF-level and full functional profiling is new**, not just faster:
  UniRef, GO, Pfam, KEGG, MetaCyc and eggNOG classification (including
  external ID cross-references, KEGG/MetaCyc compound inference, and
  pathway/module coverage) are all computed now, on top of the OGU table.
  The original's `--no-fun` meant none of this was ever produced.
- OGU-level classification itself (genome, no taxonomic ranks) is
  unchanged in scope -- both pipelines skip ranked taxonomic
  classification (phylum...species); that's not new, it was already the
  case via `--no-tax`.
- No per-sample `.sam` file is ever written -- Bowtie2 output streams
  directly into Woltka, with OGU and ORF classification running
  concurrently on that same stream.
- Work is split into parallel batches (SLURM array jobs) instead of one
  single job for the whole cohort, then automatically merged back into
  one final table set at the end -- the original never needed a merge
  step because it was never split into batches to begin with.
