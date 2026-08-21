#!/bin/bash
#SBATCH -c 16
#SBATCH --mem 128G
#SBATCH -J align_classify_batch
#SBATCH -o slurm-%x-%A_%a.out
#SBATCH -e slurm-%x-%A_%a.err
#### (submit_batches.sh overrides these two via --output/--error so logs
#### land in $LOGS_ROOT; they're only used as-is if you sbatch this file
#### directly, e.g. to retry a single failed batch.)
####
#### Streams bowtie2 output directly into woltka -- no .sam file, not even
#### a compressed one, ever touches disk. Woltka needs the alignment TWICE:
#### once for OGU (genome-level) assignment, once for ORF assignment. A
#### plain pipe can only be read once, so the combined, sample-prefixed
#### alignment stream is `tee`'d into two concurrent `woltka classify`
#### processes -- one reads from a FIFO (OGU pass), the other reads tee's
#### own stdout (ORF pass). Both consume data as bowtie2 produces it, so
#### this is still fully streamed: the FIFO is just live plumbing between
#### two processes, nothing accumulates on disk. Taxonomic rank
#### classification (phylum...species) is still deliberately excluded --
#### that needs its own separate pass too, and OGU already covers what
#### this pipeline needs. If you want ranks back, see bowtie2_batch.sh /
#### woltka_batch.sh / sam2bam_batch.sh (the .sam.gz + wol2sop.sh pipeline)
#### instead, or extend the tee below with a third branch.
####
#### How it works: each sample's bowtie2 alignment stream is piped through
#### awk, which prefixes every read ID with "<sample>_". All samples in the
#### batch are concatenated into ONE stream, which is teed to the OGU and
#### ORF classify calls. Woltka's stdin mode demultiplexes automatically on
#### the first underscore in the read ID, splitting the combined stream
#### back into one column per sample in each output table -- same
#### one-table-per-batch result as the old directory-based approach.
####
#### Below that, the functional cascade reproduces wol2sop.sh's *full*
#### feature set (everything except taxonomic rank classification): all six
#### default databases (uniref, go, pfam, kegg, metacyc, eggnog), external
#### ID cross-references (--idmaps), chemical compound inference for KEGG
#### and MetaCyc (--compound), pathway/module coverage (--coverage), and
#### gene-length normalization of the ORF table (--funnrm/--rpk). Each is
#### behind its own toggle below (IDMAPS/COMPOUND/COVERAGE/FUNNRM/TAXNRM) so
#### you can switch any of them back off without touching the cascade.
####
#### CAVEATS (read before relying on this):
#### - Sample IDs in the manifest must NOT contain "_". Demux splits on the
####   FIRST underscore; a sample ID with one in it silently misattributes
####   reads to the wrong column. This script aborts the batch if it finds
####   one -- rename the sample or fall back to the file-based pipeline.
#### - No per-sample resumability. The old bowtie2_batch.sh could skip
####   samples whose .sam already existed. Here there's nothing on disk to
####   resume from: if this job dies, the whole batch re-aligns from
####   scratch on retry.
#### - This one job now does alignment AND two concurrent classify passes
####   (connected via tee/FIFO), so it needs enough cores/memory for all
####   three at once. Tune -c/--mem above if you see contention; the
####   functional cascade below runs after the pipe finishes, not
####   concurrently with it.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
BATCH_ID=$SLURM_ARRAY_TASK_ID  # 1-based batch index passed via job array
BATCH_DIR=$(printf "Batch%02d" "$BATCH_ID")
cd "$SLURM_SUBMIT_DIR"

# OUTPUT_ROOT (from config.sh) defaults to ".", i.e. output lands under
# $SLURM_SUBMIT_DIR as before. readlink -m resolves it whether it's given
# as relative (to the submit dir) or absolute.
OUTPUT_ROOT_ABS=$(readlink -m "$OUTPUT_ROOT")
OUT="$OUTPUT_ROOT_ABS/$BATCH_DIR/Output"
mkdir -p "$OUT"

WOL2DB_ABS=$(readlink -m "$WOL2DB")
FMT=biom   # or tsv; matches wol2sop.sh's -f option
[ "$FMT" == biom ] && ALTFMT= || ALTFMT="--to-tsv"

# Toggles mirroring wol2sop.sh's optional flags. Set any of these to 1 to
# switch that piece on, or back to empty ("") to switch it off.
IDMAPS=1
COMPOUND=1
COVERAGE=1
# Off by default: unlike IDMAPS/COMPOUND/COVERAGE (which only add extra
# files), FUNNRM changes the actual values in orf.$FMT itself -- it
# normalizes read counts to reads-per-kilobase of gene length instead of
# raw counts. Your original wol2sop.sh run never passed --funnrm, so this
# is off to match those raw counts exactly. Set to 1 any time you want
# gene-length-normalized counts instead.
FUNNRM=
[ -n "$FUNNRM" ] && FUNNRM_ARGS="--sizes . --scale 1k --digits 3" || FUNNRM_ARGS=

# Genome-size normalization of the OGU table (wol2sop.sh's --taxnrm). Off
# by default: it needs $WOL2DB_ABS/genomes/length.map, and this pipeline
# deliberately skips downloading genomes/ (only needed for this one flag).
# Set to 1 and populate wol2/genomes/ if you want it.
TAXNRM=
[ -n "$TAXNRM" ] && TAXNRM_ARGS="--sizes $WOL2DB_ABS/genomes/length.map --scale 1M" || TAXNRM_ARGS=

# Reserve some of the allocated cores for the two concurrent classify
# processes (OGU + ORF) and tee, so they don't starve bowtie2 or each other.
BT2_THREADS=12

START=$(( (BATCH_ID - 1) * BATCH_SIZE + 1 ))
END=$(( START + BATCH_SIZE - 1 ))
mapfile -t ROWS < <(tail -n +2 "$MANIFEST")

# Guard: demultiplexing splits on the FIRST underscore in the read ID, so a
# sample ID containing "_" would corrupt sample assignment silently.
for ((i=START; i<=END && i<=${#ROWS[@]}; i++)); do
  SAMPLE=$(cut -f1 <<< "${ROWS[$((i - 1))]}")
  if [[ "$SAMPLE" == *_* ]]; then
    echo "ERROR: sample ID '$SAMPLE' contains '_', which breaks stdin demultiplexing. Rename it in the manifest, or use the file-based pipeline for this batch." >&2
    exit 1
  fi
done

echo "# Streaming batch $BATCH_ID ($BATCH_DIR) directly into woltka (OGU + ORF, no intermediate alignment file)."
echo "# "$(date)

##################################################
# Alignment + OGU + ORF classify (two passes,   #
# run concurrently via tee, still streamed)      #
##################################################
FIFO=$(mktemp -u "$SLURM_SUBMIT_DIR/.ogu_fifo.${BATCH_ID}.XXXXXX")
mkfifo "$FIFO"
trap 'rm -f "$FIFO"' EXIT

# Start the OGU classify process first. It blocks opening the FIFO for
# reading until tee (below) opens it for writing, so there's no race.
woltka classify \
    --input - \
    $ALTFMT \
    $TAXNRM_ARGS \
    --output "$OUT/output.$FMT" \
    < "$FIFO" &
OGU_PID=$!

# set +e around this block: we want to capture both the ORF pipeline's and
# the OGU process's exit status ourselves before deciding to abort, rather
# than have `set -e` kill the script on the first failure and leave the
# other process orphaned.
set +e
{
  for ((i=START; i<=END && i<=${#ROWS[@]}; i++)); do
    ROW="${ROWS[$((i - 1))]}"
    SAMPLE=$(cut -f1 <<< "$ROW")
    READS=$(cut -f2 <<< "$ROW")
    echo "# Aligning $SAMPLE..." >&2
    bowtie2 -x "$DB" -U "$READS" -p "$BT2_THREADS" --no-hd \
      | awk -v s="$SAMPLE" 'BEGIN{FS=OFS="\t"} {$1=s"_"$1; print}'
  done
} | tee "$FIFO" | woltka classify \
    --input - \
    $ALTFMT \
    $FUNNRM_ARGS \
    --coords "$WOL2DB_ABS/proteins/coords.txt.xz" \
    --output "$OUT/orf.$FMT"
ORF_STATUS=$?
wait "$OGU_PID"
OGU_STATUS=$?
set -e

rm -f "$FIFO"
trap - EXIT

if [[ $ORF_STATUS -ne 0 ]]; then
  echo "ERROR: alignment/ORF classification failed (status $ORF_STATUS)." >&2
  exit 1
fi
if [[ $OGU_STATUS -ne 0 ]]; then
  echo "ERROR: OGU classification failed (status $OGU_STATUS)." >&2
  exit 1
fi

echo "# OGU assignment completed (output.$FMT)."
echo "# ORF assignment completed (orf.$FMT)."
echo "# "$(date)

#######################
# Functional analysis #
#######################
# Mirrors wol2sop.sh's functional cascade (default databases: uniref, go,
# pfam, kegg, metacyc, eggnog -- "all versions"), just reusing the orf.$FMT
# table produced above instead of re-deriving it from $input.
echo "# Functional analysis started."
cd "$OUT"

##########
# UniRef #
##########
ur="$WOL2DB_ABS/function/uniref"
mkdir -p uniref
cd uniref
# UniRef entries (combined UniRef90 + UniRef50)
woltka collapse -m "$ur/orf-to-uniref.map.xz" -n "$ur/uniref_name.txt.xz" \
  -i ../orf.$FMT -o uniref.$FMT
# external databases
if [[ -n "$IDMAPS" ]]; then
  woltka collapse -m "$ur/idmaps/BioCyc.map.xz" \
    -i uniref.$FMT -o biocyc.$FMT
  woltka collapse -m "$ur/idmaps/eggNOG.map.xz" \
    -i uniref.$FMT -o eggnog.$FMT
  woltka collapse -m "$ur/idmaps/GeneID.map.xz" \
    -i uniref.$FMT -o geneid.$FMT
  woltka collapse -m "$ur/idmaps/Gene_Name.map.xz" \
    -i uniref.$FMT -o gene.$FMT
  woltka collapse -m "$ur/idmaps/OMA.map.xz" \
    -i uniref.$FMT -o oma.$FMT
  woltka collapse -m "$ur/idmaps/OrthoDB.map.xz" \
    -i uniref.$FMT -o orthodb.$FMT
  woltka collapse -m "$ur/idmaps/PATRIC.map.xz" \
    -i uniref.$FMT -o patric.$FMT
  woltka collapse -m "$ur/idmaps/RefSeq.map.xz" \
    -i uniref.$FMT -o refseq.$FMT
  woltka collapse -m "$ur/idmaps/STRING.map.xz" \
    -i uniref.$FMT -o string.$FMT
fi
cd ..

######
# GO #
######
go="$WOL2DB_ABS/function/go"
mkdir -p go
cd go
# UniRef to GO (by domain)
for domain in all component function process; do
  woltka collapse -m "$go/uniref/$domain.map.xz" -n "$go/go_name.txt" \
    -i ../uniref/uniref.$FMT -o $domain.$FMT
  # GO slim (generic)
  woltka collapse -m "$go/generic/$domain.map" -n "$go/go_name.txt" \
    -i $domain.$FMT -o $domain.generic.$FMT
done
# external databases
if [[ -n "$IDMAPS" ]]; then
  woltka collapse -m "$go/idmaps/ec.map" \
    -i all.$FMT -o ec.$FMT
  woltka collapse -m "$go/idmaps/kegg.map" \
    -i all.$FMT -o kegg.$FMT
  woltka collapse -m "$go/idmaps/metacyc.map" \
    -i all.$FMT -o metacyc.$FMT
  woltka collapse -m "$go/idmaps/reactome.map" \
    -i all.$FMT -o reactome.$FMT
  woltka collapse -m "$go/idmaps/rhea.map" \
    -i all.$FMT -o rhea.$FMT
fi
cd ..

########
# Pfam #
########
pf="$WOL2DB_ABS/function/pfam"
mkdir -p pfam
cd pfam
# ORF to Pfam
woltka collapse -m "$pf/orf-to-pfam.map.xz" -n "$pf/pfam_name.txt" \
  -i ../orf.$FMT -o pfam.$FMT
# Pfam to clan
woltka collapse -m "$pf/pfam-to-clan.map" -n "$pf/clan_name.txt" \
  -i pfam.$FMT -o clan.$FMT
# external databases
if [[ -n "$IDMAPS" ]]; then
  # Pfam to InterPro
  woltka collapse -m "$pf/pfam-to-interpro.map" \
    -i pfam.$FMT -o interpro.$FMT
  # Pfam to GO
  mkdir -p go
  for domain in all component function process; do
    woltka collapse -m "$pf/pfam-to-go/$domain.map" \
      -i pfam.$FMT -o go/$domain.$FMT
  done
fi
cd ..

########
# KEGG #
########
ke="$WOL2DB_ABS/function/kegg"
mkdir -p kegg
cd kegg
# ORF to KO
woltka collapse -m "$ke/orf-to-ko.map.xz" -n "$ke/ko_name.txt" \
  -i ../orf.$FMT -o ko.$FMT
# KO to EC
woltka collapse -m "$ke/ko-to-ec.map" \
  -i ko.$FMT -o ec.$FMT
# main cascade
# KO to reaction
woltka collapse -m "$ke/ko-to-reaction.map" -n "$ke/reaction_name.txt" \
  -i ko.$FMT -o reaction.$FMT
# reaction to module
woltka collapse -m "$ke/reaction-to-module.map" -n "$ke/module_name.txt" \
  -i reaction.$FMT -o module.$FMT
# module to pathway
woltka collapse -m "$ke/module-to-pathway.map" -n "$ke/pathway_name.txt" \
  -i module.$FMT -o pathway.$FMT
# classes
# reaction to rclass
woltka collapse -m "$ke/reaction-to-rclass.map" -n "$ke/rclass_name.txt" \
  -i reaction.$FMT -o rclass.$FMT
# module class
woltka collapse -m "$ke/module-to-class.map" \
  -i module.$FMT -o module_class.$FMT
# pathway class
woltka collapse -m "$ke/pathway-to-class.map" \
  -i pathway.$FMT -o pathway_class.$FMT
# KO to disease
woltka collapse -m "$ke/ko-to-disease.map" -n "$ke/disease_name.txt" \
  -i ko.$FMT -o disease.$FMT
# compound (incl. glycan and drug)
if [[ -n "$COMPOUND" ]]; then
  for side in left right; do
    woltka collapse -m "$ke/reaction-to-${side}_compound.map" -n "$ke/compound_name.txt" \
      -i reaction.$FMT -o ${side}_compound.$FMT
  done
  woltka merge -i left_compound.$FMT -i right_compound.$FMT -o compound.$FMT
fi
# external databases
if [[ -n "$IDMAPS" ]]; then
  # KO to GO
  woltka collapse -m "$ke/ko-to-go.map" \
    -i ko.$FMT -o go.$FMT
  # KO to COG
  woltka collapse -m "$ke/ko-to-cog.map" \
    -i ko.$FMT -o cog.$FMT
fi
# coverage
if [[ -n "$COVERAGE" ]]; then
  # module coverage by reaction
  woltka coverage -m "$ke/module-to-reaction.map" \
    -i reaction.$FMT -o module_coverage.$FMT
  # pathway coverage by module
  woltka coverage -m "$ke/pathway-to-module.map" \
    -i module.$FMT -o pathway_coverage.$FMT
fi
cd ..

###########
# MetaCyc #
###########
mc="$WOL2DB_ABS/function/metacyc"
mkdir -p metacyc
cd metacyc
# ORF to protein
woltka collapse -m "$mc/orf-to-protein.map.xz" -n "$mc/protein_name.txt" \
  -i ../orf.$FMT -o protein.$FMT
# main cascade
# protein to enzrxn (enzymatic reaction)
woltka collapse -m "$mc/protein-to-enzrxn.map" -n "$mc/enzrxn_name.txt" \
  -i protein.$FMT -o enzrxn.$FMT
# enzrxn to reaction
woltka collapse -m "$mc/enzrxn-to-reaction.map" -n "$mc/reaction_name.txt" \
  -i enzrxn.$FMT -o reaction.$FMT
# reaction to pathway
woltka collapse -m "$mc/reaction-to-pathway.map" -n "$mc/pathway_name.txt" \
  -i reaction.$FMT -o pathway.$FMT
# pathway to super pathway
woltka collapse -m "$mc/pathway-to-super_pathway.map" -n "$mc/pathway_name.txt" \
  -i pathway.$FMT -o super_pathway.$FMT
# super pathway (or pathway) to pathway type
woltka collapse -m "$mc/pathway_type.txt" -n "$mc/all_class_name.txt" \
  -i super_pathway.$FMT -o pathway_type.$FMT
# branches
# protein to gene
woltka collapse -m "$mc/protein-to-gene.map" -n "$mc/gene_name.txt" \
  -i protein.$FMT -o gene.$FMT
# enzrxn to regulation
woltka collapse -m "$mc/enzrxn-to-regulation.map" \
  -i enzrxn.$FMT -o regulation.$FMT
# regulation to regulator
woltka collapse -m "$mc/regulation-to-regulator.map" -n "$mc/compound_name.txt" \
  -i regulation.$FMT -o regulator.$FMT
# compound
if [[ -n "$COMPOUND" ]]; then
  # reaction to compound (left and right)
  for side in left right; do
    woltka collapse -m "$mc/reaction-to-${side}_compound.map" -n "$mc/compound_name.txt" \
      -i reaction.$FMT -o ${side}_compound.$FMT
    # compound type
    woltka collapse -m "$mc/compound_type.txt" -n "$mc/all_class_name.txt" \
      -i ${side}_compound.$FMT -o ${side}_compound_type.$FMT
  done
  # compound and type (both sides)
  woltka merge -i left_compound.$FMT -i right_compound.$FMT -o compound.$FMT
  woltka merge -i left_compound_type.$FMT -i right_compound_type.$FMT -o compound_type.$FMT
fi
# coverage
if [[ -n "$COVERAGE" ]]; then
  # pathway coverage (by reaction)
  woltka coverage -m "$mc/pathway-to-reaction_list.map" \
    -i reaction.$FMT -o pathway_coverage.$FMT
fi
# external databases
if [[ -n "$IDMAPS" ]]; then
  # protein to go
  woltka collapse -m "$mc/protein-to-go.map" \
    -i protein.$FMT -o go.$FMT
  # reaction to EC
  woltka collapse -m "$mc/reaction-to-ec.map" \
    -i reaction.$FMT -o ec.$FMT
fi
cd ..

##########
# eggNOG #
##########
en="$WOL2DB_ABS/function/eggnog"
mkdir -p eggnog
cd eggnog
# ORF to seed orthologs
woltka collapse -m "$en/orf-to-seed.map.xz" \
  -i ../orf.$FMT -o seed.$FMT
# seed to gene
woltka collapse -m "$en/seed-to-gene.map.xz" \
  -i seed.$FMT -o gene.$FMT
# seed to orthologous group (OG)
woltka collapse -m "$en/seed-to-og.map.xz" -n "$en/og_description.txt.xz" \
  -i seed.$FMT -o og.$FMT
# OG to category
woltka collapse -m "$en/og-to-category.map.xz" -n "$en/cog_category.txt" \
  -i og.$FMT -o category.$FMT
# external databases
if [[ -n "$IDMAPS" ]]; then
  # KEGG catalogs
  mkdir -p kegg
  woltka collapse -m "$en/idmaps/KEGG_ko.map.xz" \
    -i seed.$FMT -o kegg/ko.$FMT
  woltka collapse -m "$en/idmaps/KEGG_Pathway.map.xz" \
    -i seed.$FMT -o kegg/pathway.$FMT
  woltka collapse -m "$en/idmaps/KEGG_Module.map.xz" \
    -i seed.$FMT -o kegg/module.$FMT
  woltka collapse -m "$en/idmaps/KEGG_Reaction.map.xz" \
    -i seed.$FMT -o kegg/reaction.$FMT
  woltka collapse -m "$en/idmaps/KEGG_rclass.map.xz" \
    -i seed.$FMT -o kegg/rclass.$FMT
  woltka collapse -m "$en/idmaps/BRITE.map.xz" \
    -i seed.$FMT -o kegg/brite.$FMT
  woltka collapse -m "$en/idmaps/KEGG_TC.map.xz" \
    -i seed.$FMT -o kegg/tc.$FMT
  # other databases
  woltka collapse -m "$en/idmaps/GOs.map.xz" \
    -i seed.$FMT -o go.$FMT
  woltka collapse -m "$en/idmaps/EC.map.xz" \
    -i seed.$FMT -o ec.$FMT
  woltka collapse -m "$en/idmaps/CAZy.map.xz" \
    -i seed.$FMT -o cazy.$FMT
  woltka collapse -m "$en/idmaps/BiGG_Reaction.map.xz" \
    -i seed.$FMT -o bigg.$FMT
  woltka collapse -m "$en/idmaps/PFAMs.map.xz" \
    -i seed.$FMT -o pfam.$FMT
fi
cd ..

echo "# Functional analysis completed."
echo "# "$(date)
