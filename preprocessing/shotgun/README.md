# Shotgun metagenomics preprocessing

This directory contains two preprocessing pipelines for shotgun metagenomic
sequencing data, both designed to run on an HPC cluster using SLURM.

## `PreProcessingShotgun/`

The original pipeline used to generate the shotgun metagenomics results
reported in the manuscript (Kuiper et al. 2026). Runs as a single,
non-batched SLURM job for the whole cohort: every sample is aligned with
Bowtie2 in one job, each written to disk as a `.sam` file, then classified
once with `wol2sop.sh -d ./wol2 -i align -o Output -f biom --no-tax
--no-fun` -- OGU-level (genome) classification only; `--no-fun` means no
ORF assignment or functional profiling is done at all. Kept as-is for
reproducibility of the published results -- please don't modify it if you
need to reproduce the manuscript's numbers exactly.

## `PreProcessingShotgunUpdate/`

A faster, more disk/memory-efficient version of the same pipeline, for
processing additional cohort data going forward -- and one that adds ORF
and functional profiling, which the original never computed. Main
differences from `PreProcessingShotgun/`:

- **ORF-level and full functional profiling is new**, not just faster:
  UniRef, GO, Pfam, KEGG, MetaCyc and eggNOG classification (including
  external ID cross-references, KEGG/MetaCyc compound inference, and
  pathway/module coverage) are all computed now. The original ran with
  `--no-fun` and only ever produced the OGU table.
- OGU-level classification itself (genome, no taxonomic ranks) is
  unchanged in scope from the original -- both skip ranked taxonomic
  classification (phylum...species); that was already the case via
  `--no-tax`.
- No per-sample `.sam` file is ever written -- Bowtie2 output streams
  directly into Woltka, with OGU and ORF classification running
  concurrently on that same stream (fanned out with `tee`).
- Work is split into parallel batches (SLURM array jobs) instead of one
  single job for the whole cohort, then automatically merged back into
  one final table set at the end -- the original never needed a merge
  step because it was never split into batches to begin with.

See `PreProcessingShotgunUpdate/README.md` for the full setup and usage
walkthrough.

## Which one should I use?

- **Reproducing the published manuscript results exactly:** use
  `PreProcessingShotgun/`.
- **Processing new or additional shotgun samples going forward:** use
  `PreProcessingShotgunUpdate/`.

## After preprocessing: GreenGenes2 harmonization

Whichever pipeline you use, its output (a `woltka.biom.qza`) still needs
to be called against GreenGenes2 -- this is a separate, shared follow-up
step that lives here in `shotgun/` (not inside either preprocessing
subfolder), since it's identical regardless of which pipeline produced
the `.qza`. `woltka.biom.qza` should sit in this folder (or wherever you
run these two scripts from) before you start -- see
`PreProcessingShotgunUpdate/README.md` for how that pipeline gets it there
automatically. Then run, **in this order**:

1. [`closed_referenceShotgun.sbatch`](closed_referenceShotgun.sbatch) --
   calls `woltka.biom.qza` against the GreenGenes2 backbone (`qiime
   greengenes2 filter-features`), producing `feature.table.gg2-2022.10.qza`.
2. [`taxonomic_table_Shotgun.sbatch`](taxonomic_table_Shotgun.sbatch) --
   derives the corresponding taxonomy (`qiime greengenes2
   taxonomy-from-table`) from that filtered table, producing
   `df.gg2.taxonomy.qza`.

These currently call GreenGenes2 release **2022.10** (`--i-reference
2022.10.taxonomy.asv.nwk.qza`, matching the output filenames above). To
use a newer release (e.g. 2024.09) instead, edit the reference file and
output filenames in both `closed_referenceShotgun.sbatch` and
`taxonomic_table_Shotgun.sbatch` accordingly.

See the top-level repository README's "GreenGenes2" section for full
details. For shotgun data, the only prerequisite is the
`2022.10.taxonomy.asv.nwk.qza` reference already downloaded during the
initial QIIME2 setup at the top of that README (the backbone/`get_repset.py`
prerequisites there apply only to non-V4 16S data, not shotgun).

## Before either pipeline: quality control

Both pipelines expect already quality-controlled (kneaddata-processed)
reads as input -- QC is expected to be done separately, before either
pipeline.
