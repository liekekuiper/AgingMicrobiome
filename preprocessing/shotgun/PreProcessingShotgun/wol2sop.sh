#!/bin/bash

#######################################
#
# Standard workflow for WoL2 + Woltka
#
#######################################

# Author:  Qiyun Zhu
# License: BSD-3-Clause
# Version: 0.0.1-dev
# Email: qiyunzhu@gmail.com

# Last updated: 2022-12-25

# This workflow requires Woltka v0.1.5 or above.
# See: https://github.com/qiyunzhu/woltka

########
# Help #
########

function help() {
  cat << HELP

  Standard workflow for WoL2 + Woltka

  Usage: $(basename $0) -d WOL2DB -i INPUT [OPTIONS]

  Required parameters:
    -d, --db             path to WoL2 database directory
    -i, --input          path to input alignment file or directory

  Input and output:
    -o, --output         path to output directory (default: current)
    -f, --outfmt         output table format: biom (default) or tsv
    -e, --filext         alignment filename extension

  Taxonomic analysis:
    --no-tax             skip taxonomic analysis (only infer OGUs)
    --taxlvs             taxonomic ranks to classify to
                         (default: phylum,class,order,family,genus,species)
    --taxfmt             taxonomy tree format (taxdump or lineage)
    --taxnrm             normalize read count by genome size (Mbp)

  Functional analysis:
    --no-fun             skip functional analysis (only infer ORFs)
    --fundbs             function databases to use
                         (default: uniref,go,pfam,kegg,metacyc,eggnog)
    --funnrm, --rpk      normalize read count by gene size (kbp); result is
                         RPK (reads per kilobase)
    --compound           infer chemical compounds (for KEGG & MetaCyc)
    --coverage           calculate pathway coverage
    --idmaps             generate mappings to external databases

  Help information:
    -h, --help           show this help message and exit

HELP
}


##############
# Parameters #
##############

# exit on error

set -e

# show help info

if [[ $# -eq 0 ]]; then
  help
  exit 1
fi

echo "# Standard workflow for WoL2 + Woltka"
echo "# "$(date)

# read paramters

while [[ $# -gt 0 ]]; do
  case $1 in
    -h | --help  ) help; exit 1 ;;

    -d | --db    ) db="$2";     shift ;;
    -i | --input ) input="$2";  shift ;;

    -o | --output) output="$2"; shift ;;
    -f | --outfmt) fmt="$2";    shift ;;
    -e | --filext) ext="$2";    shift ;;

    --no-tax) notax=1                 ;;
    --taxlvs) taxlvs="$2";      shift ;;
    --taxfmt) taxfmt="$2";      shift ;;
    --taxnrm) taxnrm=1                ;;

    --no-fun) nofun=1                 ;;
    --fundbs) fundbs="$2";      shift ;;
    --funnrm | --rpk) funnrm=1        ;;
    --compound) compound=1            ;;
    --coverage) coverage=1            ;;
    --idmaps) idmaps=1                ;;

    *) echo "Invalid parameter: $1."; exit 1 ;;
  esac
  shift
done


##############
# Validation #
##############

# woltka executable

if ! command -v woltka &> /dev/null; then
  echo "Missing Woltka executable."
  exit
fi

ver=$(woltka --version | cut -f3- -d' ')
echo "# Woltka version: "$ver

# database directory

if [[ -z "$db" ]]; then
  echo "Missing parameter: -d|--db"
  exit 1
fi

if [[ ! -d "$db" ]]; then
  echo "Database directory $db does not exist."
  exit 1
fi

db=$(readlink -m $db)  # get absolute path

echo "# Database directory: $db"

# input file / directory

if [[ -z "$input" ]]; then
  echo "Missing parameter: -i|--input"
  exit 1
fi

if [[ ! -f "$input" ]] && [[ ! -d "$input" ]]; then
  echo "Input file or directory $input does not exist."
  exit 1
fi

input=$(readlink -m $input)

if [[ -f $input ]]; then
  echo "# Input file: $input"
else
  echo "# Input directory: $input"
fi

# output directory

if [[ -z "$output" ]]; then
  output=.
elif [[ -f "$output" ]]; then
  echo "Output should be a directory not a file."
  exit 1
elif [[ ! -d "$output" ]]; then
  mkdir -p $output
fi

output=$(readlink -m $output)

echo "# Output directory: $output"

# other parameters

[ -z "$fmt" ] && fmt=biom

[ "$fmt" == biom ] && altfmt= || altfmt="--to-tsv"

[ -z "$ext" ] && filext= || filext="--filext $ext"

[ -z "$taxfmt" ] && taxfmt=taxdump

[ -z "$taxlvs" ] && taxlvs=phylum,class,order,family,genus,species

[ -z "$fundbs" ] && fundbs=uniref,go,pfam,kegg,metacyc,eggnog

[ -z "$taxnrm" ] && taxnrm= || taxnrm="--sizes $db/genomes/length.map --scale 1M"

[ -z "$funnrm" ] && funnrm= || funnrm="--sizes . --scale 1k --digits 3"


#######
# SOP #
#######

cwd=$(pwd)
cd $output

# taxonomic classification

if [[ -z "$notax" ]]; then
  echo "# Taxonomic analysis started."
  echo "# Will classify to OGU (none), free-rank, and the following ranks:"
  echo "#   "$taxlvs | sed 's/,/, /g'

  if [ "$taxfmt" == taxdump ]; then
    woltka classify $filext $altfmt $taxnrm \
      --input  $input \
      --map    $db/taxonomy/taxid.map \
      --nodes  $db/taxonomy/nodes.dmp \
      --names  $db/taxonomy/names.dmp \
      --rank   none,free,$taxlvs \
      --output .

  else
    woltka classify $filext $altfmt $taxnrm \
      --input   $input \
      --lineage $db/taxonomy/lineages.txt \
      --rank    none,free,$taxlvs \
      --output .
  fi

  echo "# Taxonomic analysis completed."

# OGU only

else
  echo "# OGU assignment started."

  woltka classify $filext $altfmt $taxnrm \
    --input $input \
    --output output.biom

  echo "# OGU assignment completed."
fi

#mv none.$fmt ogu.$fmt


###################
# Read/gene match #
###################

echo "# ORF assignment started."

woltka classify $filext $funnrm \
  --input  $input \
  --coords $db/proteins/coords.txt.xz \
  --output orf.$fmt

echo "# ORF assignment completed."

# skip functional analysis

if [[ ! -z "$nofun" ]]; then
  cd $cwd
  echo $(date)
  exit 0
fi

echo "# Functional analysis started."
echo "# Will classify using the following databases:"
echo "#   "$fundbs | sed 's/,/, /g'

fundbs=(${fundbs//,/ })  # convert to array


##########
# UniRef #
##########

if [[ " ${fundbs[*]} " =~ " uniref " ]] || [[ " ${fundbs[*]} " =~ " go " ]]; then
  echo "# Classifying using UniRef..."

  ur=$db/function/uniref
  mkdir -p uniref
  cd uniref

  # UniRef entries (combined UniRef90 + UniRef50)
  woltka collapse -m $ur/orf-to-uniref.map.xz -n $ur/uniref_name.txt.xz \
    -i ../orf.$fmt -o uniref.$fmt

  # external databases
  if [[ ! -z "$idmaps" ]] && [[ " ${fundbs[*]} " =~ " uniref " ]]; then

    woltka collapse -m $ur/idmaps/BioCyc.map.xz \
      -i uniref.$fmt -o biocyc.$fmt

    woltka collapse -m $ur/idmaps/eggNOG.map.xz \
      -i uniref.$fmt -o eggnog.$fmt

    woltka collapse -m $ur/idmaps/GeneID.map.xz \
      -i uniref.$fmt -o geneid.$fmt

    woltka collapse -m $ur/idmaps/Gene_Name.map.xz \
      -i uniref.$fmt -o gene.$fmt

    woltka collapse -m $ur/idmaps/OMA.map.xz \
      -i uniref.$fmt -o oma.$fmt

    woltka collapse -m $ur/idmaps/OrthoDB.map.xz \
      -i uniref.$fmt -o orthodb.$fmt

    woltka collapse -m $ur/idmaps/PATRIC.map.xz \
      -i uniref.$fmt -o patric.$fmt

    woltka collapse -m $ur/idmaps/RefSeq.map.xz \
      -i uniref.$fmt -o refseq.$fmt

    woltka collapse -m $ur/idmaps/STRING.map.xz \
      -i uniref.$fmt -o string.$fmt

  fi
  cd ..
fi


######
# GO #
######


if [[ " ${fundbs[*]} " =~ " go " ]]; then
  echo "# Classifying using GO..."

  go=$db/function/go
  mkdir -p go
  cd go

  # UniRef to GO (by domain)
  for domain in all component function process; do
    woltka collapse -m $go/uniref/$domain.map.xz -n $go/go_name.txt \
      -i ../uniref/uniref.$fmt -o $domain.$fmt

    # GO slim (generic)
    woltka collapse -m $go/generic/$domain.map -n $go/go_name.txt \
      -i $domain.$fmt -o $domain.generic.$fmt
  done

  # external databases
  if [[ ! -z "$idmaps" ]]; then

    woltka collapse -m $go/idmaps/ec.map \
      -i all.$fmt -o ec.$fmt

    woltka collapse -m $go/idmaps/kegg.map \
      -i all.$fmt -o kegg.$fmt

    woltka collapse -m $go/idmaps/metacyc.map \
      -i all.$fmt -o metacyc.$fmt

    woltka collapse -m $go/idmaps/reactome.map \
      -i all.$fmt -o reactome.$fmt

    woltka collapse -m $go/idmaps/rhea.map \
      -i all.$fmt -o rhea.$fmt

  fi
  cd ..
fi


########
# Pfam #
########

if [[ " ${fundbs[*]} " =~ " pfam " ]]; then
  echo "# Classifying using Pfam..."

  pf=$db/function/pfam
  mkdir -p pfam
  cd pfam

  # ORF to Pfam
  woltka collapse -m $pf/orf-to-pfam.map.xz -n $pf/pfam_name.txt \
    -i ../orf.$fmt -o pfam.$fmt

  # Pfam to clan
  woltka collapse -m $pf/pfam-to-clan.map -n $pf/clan_name.txt \
    -i pfam.$fmt -o clan.$fmt

  # external databases
  if [[ ! -z "$idmaps" ]]; then

    # Pfam to InterPro
    woltka collapse -m $pf/pfam-to-interpro.map \
      -i pfam.$fmt -o interpro.$fmt

    # Pfam to GO
    mkdir -p go
    for domain in all component function process; do
      woltka collapse -m $pf/pfam-to-go/$domain.map \
        -i pfam.$fmt -o go/$domain.$fmt
    done

  fi
  cd ..
fi


########
# KEGG #
########

if [[ " ${fundbs[*]} " =~ " kegg " ]]; then
  echo "# Classifying using KEGG..."

  ke=$db/function/kegg
  mkdir -p kegg
  cd kegg

  # ORF to KO
  woltka collapse -m $ke/orf-to-ko.map.xz -n $ke/ko_name.txt \
    -i ../orf.$fmt -o ko.$fmt

  # KO to EC
  woltka collapse -m $ke/ko-to-ec.map \
    -i ko.$fmt -o ec.$fmt

  # main cascade
  # KO to reaction
  woltka collapse -m $ke/ko-to-reaction.map -n $ke/reaction_name.txt \
    -i ko.$fmt -o reaction.$fmt

  # reaction to module
  woltka collapse -m $ke/reaction-to-module.map -n $ke/module_name.txt \
    -i reaction.$fmt -o module.$fmt

  # module to pathway
  woltka collapse -m $ke/module-to-pathway.map -n $ke/pathway_name.txt \
    -i module.$fmt -o pathway.$fmt

  # classes
  # reaction to rclass
  woltka collapse -m $ke/reaction-to-rclass.map -n $ke/rclass_name.txt \
    -i reaction.$fmt -o rclass.$fmt

  # module class
  woltka collapse -m $ke/module-to-class.map \
    -i module.$fmt -o module_class.$fmt

  # pathway class
  woltka collapse -m $ke/pathway-to-class.map \
    -i pathway.$fmt -o pathway_class.$fmt

  # KO to disease
  woltka collapse -m $ke/ko-to-disease.map -n $ke/disease_name.txt \
    -i ko.$fmt -o disease.$fmt

  # compound (incl. glycan and drug)
  if [[ ! -z "$compound" ]]; then

    for side in left right; do
      woltka collapse -m $ke/reaction-to-${side}_compound.map -n $ke/compound_name.txt \
        -i reaction.$fmt -o ${side}_compound.$fmt
    done

    woltka merge -i left_compound.$fmt -i right_compound.$fmt -o compound.$fmt

  fi

  # external databases
  if [[ ! -z "$idmaps" ]]; then

    # KO to GO
    woltka collapse -m $ke/ko-to-go.map \
      -i ko.$fmt -o go.$fmt

    # KO to COG
    woltka collapse -m $ke/ko-to-cog.map \
      -i ko.$fmt -o cog.$fmt

  fi

  # coverage
  if [[ ! -z "$coverage" ]]; then

    # module coverage by reaction
    woltka coverage -m $ke/module-to-reaction.map \
      -i reaction.$fmt -o module_coverage.$fmt

    # pathway coverage by module
    woltka coverage -m $ke/pathway-to-module.map \
      -i module.$fmt -o pathway_coverage.$fmt

  fi
  cd ..
fi


###########
# MetaCyc #
###########

if [[ " ${fundbs[*]} " =~ " metacyc " ]]; then
  echo "# Classifying using MetaCyc..."

  mc=$db/function/metacyc
  mkdir -p metacyc
  cd metacyc

  # ORF to protein
  woltka collapse -m $mc/orf-to-protein.map.xz -n $mc/protein_name.txt \
    -i ../orf.$fmt -o protein.$fmt

  # main cascade
  # protein to enzrxn (enzymatic reaction)
  woltka collapse -m $mc/protein-to-enzrxn.map -n $mc/enzrxn_name.txt \
    -i protein.$fmt -o enzrxn.$fmt

  # enzrxn to reaction
  woltka collapse -m $mc/enzrxn-to-reaction.map -n $mc/reaction_name.txt \
    -i enzrxn.$fmt -o reaction.$fmt

  # reaction to pathway
  woltka collapse -m $mc/reaction-to-pathway.map -n $mc/pathway_name.txt \
    -i reaction.$fmt -o pathway.$fmt

  # pathway to super pathway
  woltka collapse -m $mc/pathway-to-super_pathway.map -n $mc/pathway_name.txt \
    -i pathway.$fmt -o super_pathway.$fmt

  # super pathway (or pathway) to pathway type
  woltka collapse -m $mc/pathway_type.txt -n $mc/all_class_name.txt \
  -i super_pathway.$fmt -o pathway_type.$fmt

  # branches
  # protein to gene
  woltka collapse -m $mc/protein-to-gene.map -n $mc/gene_name.txt \
    -i protein.$fmt -o gene.$fmt

  # enzrxn to regulation
  woltka collapse -m $mc/enzrxn-to-regulation.map \
    -i enzrxn.$fmt -o regulation.$fmt

  # regulation to regulator
  woltka collapse -m $mc/regulation-to-regulator.map -n $mc/compound_name.txt \
    -i regulation.$fmt -o regulator.$fmt

  # compound
  if [[ ! -z "$compound" ]]; then

    # reaction to compound (left and right)
    for side in left right; do
      woltka collapse -m $mc/reaction-to-${side}_compound.map -n $mc/compound_name.txt \
        -i reaction.$fmt -o ${side}_compound.$fmt

      # compound type
      woltka collapse -m $mc/compound_type.txt -n $mc/all_class_name.txt \
        -i ${side}_compound.$fmt -o ${side}_compound_type.$fmt
    done

    # compound and type (both sides)
    woltka merge -i left_compound.$fmt -i right_compound.$fmt -o compound.$fmt
    woltka merge -i left_compound_type.$fmt -i right_compound_type.$fmt -o compound_type.$fmt

  fi

  # coverage
  if [[ ! -z "$coverage" ]]; then

    # pathway coverage (by reaction)
    woltka coverage -m $mc/pathway-to-reaction_list.map \
      -i reaction.$fmt -o pathway_coverage.$fmt

  fi

  # external databases
  if [[ ! -z "$idmaps" ]]; then

    # protein to go
    woltka collapse -m $mc/protein-to-go.map \
      -i protein.$fmt -o go.$fmt

    # reaction to EC
    woltka collapse -m $mc/reaction-to-ec.map \
      -i reaction.$fmt -o ec.$fmt

  fi
  cd ..
fi


##########
# eggNOG #
##########

if [[ " ${fundbs[*]} " =~ " eggnog " ]]; then
  echo "# Classifying using eggNOG..."

  en=$db/function/eggnog
  mkdir -p eggnog
  cd eggnog

  # ORF to seed orthologs
  woltka collapse -m $en/orf-to-seed.map.xz \
    -i ../orf.$fmt -o seed.$fmt

  # seed to gene
  woltka collapse -m $en/seed-to-gene.map.xz \
    -i seed.$fmt -o gene.$fmt

  # seed to orthologous group (OG)
  woltka collapse -m $en/seed-to-og.map.xz -n $en/og_description.txt.xz \
    -i seed.$fmt -o og.$fmt

  # OG to category
  woltka collapse -m $en/og-to-category.map.xz -n $en/cog_category.txt \
    -i og.$fmt -o category.$fmt

  # external databases
  if [[ ! -z "$idmaps" ]]; then

    # KEGG catalogs
    mkdir -p kegg

    woltka collapse -m $en/idmaps/KEGG_ko.map.xz \
      -i seed.$fmt -o kegg/ko.$fmt

    woltka collapse -m $en/idmaps/KEGG_Pathway.map.xz \
      -i seed.$fmt -o kegg/pathway.$fmt

    woltka collapse -m $en/idmaps/KEGG_Module.map.xz \
      -i seed.$fmt -o kegg/module.$fmt

    woltka collapse -m $en/idmaps/KEGG_Reaction.map.xz \
      -i seed.$fmt -o kegg/reaction.$fmt

    woltka collapse -m $en/idmaps/KEGG_rclass.map.xz \
      -i seed.$fmt -o kegg/rclass.$fmt

    woltka collapse -m $en/idmaps/BRITE.map.xz \
      -i seed.$fmt -o kegg/brite.$fmt

    woltka collapse -m $en/idmaps/KEGG_TC.map.xz \
      -i seed.$fmt -o kegg/tc.$fmt

    # other databases
    woltka collapse -m $en/idmaps/GOs.map.xz \
      -i seed.$fmt -o go.$fmt

    woltka collapse -m $en/idmaps/EC.map.xz \
      -i seed.$fmt -o ec.$fmt

    woltka collapse -m $en/idmaps/CAZy.map.xz \
      -i seed.$fmt -o cazy.$fmt

    woltka collapse -m $en/idmaps/BiGG_Reaction.map.xz \
      -i seed.$fmt -o bigg.$fmt

    woltka collapse -m $en/idmaps/PFAMs.map.xz \
      -i seed.$fmt -o pfam.$fmt

  fi
  cd ..
fi

echo "# Functional analysis completed."
cd $cwd
echo "# "$(date)

