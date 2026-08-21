# Shared configuration for all batch scripts
MANIFEST=manifest    # Path to the sample manifest file
BATCH_SIZE=10        # Number of samples per batch
DB=./db              # Path to bowtie2 database index
WOL2DB=./wol2        # Path to WoL2/Woltka database directory (contains proteins/, function/, taxonomy/)

# Where BatchNN/ and Merged/ output folders get created. Recommended:
# woltka_files/ (keeps the (many, large) per-batch and merged tables out
# of Scripts/ itself). Relative paths are resolved relative to the submit
# directory; an absolute path (e.g. /scratch/user/Microbiome/woltka_output)
# works too.
OUTPUT_ROOT=woltka_files

# Where SLURM stdout/stderr logs from all batch scripts get written.
# Created automatically by submit_batches.sh if it doesn't exist yet.
# make_bowtie2_db.sh (run manually, outside submit_batches.sh) needs this
# directory to already exist before you submit it -- see its header.
LOGS_ROOT=Logs

# Off by default: orf.biom is the raw, uncollapsed per-ORF table -- at 1347
# samples that was already 3.5M+ features, by far the heaviest thing
# merge_batches.sh has to merge (it's what OOM-killed the 32G run). All
# per-batch BatchNN/Output/orf.biom files still exist either way, this only
# controls whether a single merged orf.biom is also built. Set to 1 if you
# do need one merged ORF-level table and have the memory for it (see the
# --mem note in merge_batches.sh).
MERGE_ORF=

# On by default: after merge_batches.sh finishes, also import the merged
# output.biom into QIIME2 as a .qza (make_q2_import.sh). Set to empty ("")
# if you don't want a .qza -- submit_batches.sh/run_pipeline.sh will then
# skip that step entirely (and you won't need QIIME2_ENV below at all).
RUN_Q2_IMPORT=1

# Conda/mamba environment with QIIME2 installed, used by make_q2_import.sh
# to import the final merged output.biom as a .qza. Matches PreProcessingShotgun's
# original environment name; change if yours differs. Unused if
# RUN_Q2_IMPORT is off.
QIIME2_ENV=qiime2-2023.7
