# Aging Microbiome — Code Repository

This repository contains all analysis code accompanying the manuscript:

> **Robust stool microbiome signatures of human aging across cohorts and sequencing platforms**  
> Kuiper et al. (2026)

## Overview

This project investigates the relationship between the stool microbiome and human aging across 
multiple population-based cohorts and sequencing platforms (16S rRNA gene amplicon sequencing 
and shotgun metagenomics). We identify robust stool microbiome signatures associated with 
chronological age and the frailty index across six US and Dutch large-scale cohorts ([Doetinchem Cohort Study](https://doi.org/10.1093/ije/dym292),
[Framingham Heart Study](https://doi.org/10.1093/ije/dyv337), [Hispanic Community Health Study/Study of Latinos](https://doi.org/10.1016/j.annepidem.2010.03.015),
[Lifelines](https://doi.org/10.1093/ije/dyab257), [the Osteoporotic Fractures in Men (MrOS) Cohort](https://doi.org/10.1002/jbmr.4518), and
[the Rotterdam Study](https://doi.org/10.1007/s10654-023-01094-1)), and validate 
these findings for chronological age and mortality in an external cohort [FINRISK](https://doi.org/10.1093/ije/dyx239). 
Furthermore, we used metabolic simulations in [the Study of Health in Pomerania](https://doi.org/10.1093/ije/dyac034) to link discovered microbial composition to metabolic function.

---

## Repository structure

### `preprocessing/`
Preprocessing pipelines for raw microbiome sequencing data into [Greengenes2](https://doi.org/10.1038/s41587-023-01845-1) harmonized data,
designed to run on an HPC cluster using SLURM.
> **Note for consortium analysts:** American cohorts may use [Qiita](https://qiita.ucsd.edu/) 
  > for preprocessing. European cohorts cannot use Qiita due to GDPR restrictions on uploading 
  > participant data to external servers and should use this pipeline instead.

- `16S/` — Pipeline for 16S rRNA gene amplicon sequencing data.  
- `shotgun/` — Pipeline for shotgun metagenomics data
  - `PreProcessingUpdate/` — Scripts to process shotgun data in parallel batches using Bowtie2 
    for alignment and Woltka for taxonomic profiling, substantially reducing processing time 
    compared to running samples sequentially

---

### `downstreamanalyses/`
Cohort-level association analyses and cross-cohort harmonization checks. These scripts were 
run by analysts within each participating cohort.

- `analysiscode.py` — Main analysis script for HPC/SLURM environments
- `analysiscode_nonslurm.py` — Version for local/non-HPC use
- `microbiome_utils.py` — Shared utility functions
- `submit.sbatch` — SLURM submission scripts
- `changenames.sbatch` — Script to modify cohort-specific naming
- `cohort.txt` — Cohort name configuration
- `mods_agingmicrobiome.txt` — List of models run
- `outs_agingmicrobiome.txt` — Expected output files
- `subsets_agingmicrobiome.txt` — Subgroup definitions (see below)

**Subfolders:**
- `harmonization_allcohorts/` — Scripts to check harmonization across cohorts and platforms:
  - PCA of stool microbiome composition across cohorts
  - UpSet plots showing overlap in detected associations
- `multiplemethods/` — Comparison of two sequencing methods (16S rRNA gene amplicon sequencing 
  and shotgun metagenomics) in cohorts where both were available from the same stool samples 
  (FHS, HCHS/SOL, and MrOS):
  - PCoA-based comparison across methods
  - Venn diagrams of overlapping associations across methods

#### Subgroup definitions
Analyses were run in the following subgroups:

| Subset | Description |
|--------|-------------|
| `all` | All participants |
| `men` | Men only |
| `women` | Women only |
| `age_1` | Chronological age ≥ 18 and < 40 |
| `age_2` | Chronological age ≥ 40 and < 50 |
| `age_3` | Chronological age ≥ 50 and < 60 |
| `age_4` | Chronological age ≥ 60 and < 70 |
| `age_5` | Chronological age ≥ 70 |

---

### `metaanalyses/`
Random effects meta-analysis of cohort-specific results and external validation, run centrally 
after collecting results from all participating discovery cohorts.

- `MicrobiomeAging_functions.R` — Core R functions shared across all meta-analysis scripts
- `MetaAnalysis_Genus.R` — Random effects meta-analysis of chronological age and frailty 
  associations at the genus level
- `MetaAnalysis_Species.R` — Random effects meta-analysis of chronological age and frailty 
  associations at the species level
- `Table2_3_Fig3.R` — Code to reproduce main manuscript tables (2, 3) and Figure 3
- `AGORA2_TableFig.R` — Metabolic modeling analyses using the AGORA2 resource to 
  functionally interpret stool microbiome findings
- `externalvalidation/` — Validation of meta-analysis findings for chronological age and 
  mortality in independent external cohort FINRISK on 16S and metagenomics dataset:
  - `analysiscode_rfile.R` — Analysis script for external validation cohort
  - `microbiome_utils.R` — Utility functions for validation analyses
  - `Biomarker_TSE_phyloseq.R` — Biomarker analyses for TSE and phyloseq data
  - `ValidationResults.R` — Summary and visualization of validation results

---

### `microbiomebiom/`
An R package that provides a unified interface to compute genus-level stool microbiome risk 
scores (genusFI) from 16S rRNA gene amplicon or shotgun metagenomics data. It supports QIIME2 
feature tables (`.qza`), phyloseq objects, and TreeSummarizedExperiment objects. The package 
includes the default GenusFI beta coefficients (`betas_default`) derived from this study, but 
also accepts custom coefficients, allowing researchers to apply or validate stool microbiome 
aging scores in their own cohorts.

**Install:**
```r
# install.packages("devtools")
devtools::install_github("liekekuiper/agingmicrobiome/microbiomebiom")
```

**Quick start:**
```r
library(microbiomebiom)

# Compute genusFI score using default coefficients
res <- compute_risk_score(ps, input_type = "phyloseq")

# Or supply your own beta coefficients
my_betas <- read.csv("my_betas.csv")
res <- compute_risk_score(ps, input_type = "phyloseq", betas = my_betas)
```

See the package README in the subfolder for full documentation.

---

## Citation
If you use this code, please cite this repository.

---

## Contact

For questions about this repository or the study:  
**Lieke Kuiper** — l.m.kuiper[at]erasmusmc.nl  
Department of Internal Medicine, Erasmus MC, Rotterdam, the Netherlands
