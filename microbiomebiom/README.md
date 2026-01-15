# microbiomebiom R Package

**Compute genus-level microbiome risk scores using default or custom coefficients.**

---

## Description

`microbiomebiom` provides a unified interface to compute genusFI or other genus-based microbiome risk scores from 16S or WGS data. It supports:

* QIIME2 feature tables (`.qza`)
* phyloseq objects (`.rds` or in-memory)
* TreeSummarizedExperiment objects (`.rds` or in-memory)

The function does **not perform cohort-based selection**; users are expected to apply inclusion/exclusion criteria prior to scoring. Metadata is used only to align sample identifiers.

---

## Installation

```r
# Install devtools if not already installed
install.packages("devtools")

# Install from your GitHub repository
devtools::install_github("liekekuiper/agingmicrobiome/microbiomebiom")
```

---

## Usage

```r
library(microbiomebiom)
```

### Default betas (genusFI)

```r
res <- compute_risk_score(ps, input_type = "phyloseq")
head(res)
```

**Output:** `data.frame` with:

| sampleid | genusFI |
| -------- | ------- |
| S1       | 0.12    |
| S2       | -0.05   |

* Uses **default GenusFI beta coefficients** (`betas_default`).

### Custom betas (risk_score)

```r
my_betas <- read.csv("my_betas.csv", stringsAsFactors = FALSE)
res <- compute_risk_score(ps, input_type = "phyloseq", betas = my_betas)
head(res)
```

**Output:** `data.frame` with:

| sampleid | risk_score |
| -------- | ---------- |
| S1       | 0.24       |
| S2       | -0.10      |

**Beta CSV format:**

| RISK               | B     |
| ------------------ | ----- |
| g__Lactobacillus   | 0.15  |
| g__Bifidobacterium | -0.07 |

* `RISK` = genus name with `g__` prefix
* `B` = numeric coefficient
* Only genera present in both data and CSV are used

### Optional metadata

```r
res <- compute_risk_score(ps, input_type = "phyloseq",
                          metadata = metadata_df,
                          sampleid_col = "my_sample_id")
```

* Metadata is **only for aligning sample IDs**. If omitted, all microbiome samples are scored.

---

## Inputs

* `input`: microbiome data (QIIME2 `.qza`, phyloseq, or TSE object)
* `taxonomy_qza`: required if `input_type = "qza"`
* `metadata`: optional `data.frame` with sample IDs
* `sampleid_col`: column name in metadata with sample IDs (default: `sampleid`)
* `input_type`: one of `c("qza", "phyloseq", "tse")`
* `betas`: data.frame with columns `RISK` and `B` (default: `betas_default`)

---

## Output

* `data.frame` with columns depending on betas used:

  * Default betas: `sampleid | genusFI`
  * Custom betas: `sampleid | risk_score`

Additionally, it returns all samples for which microbiome data are available or aligned to metadata.

---

## Example Workflow

```r
library(microbiomebiom)
library(phyloseq)

# Load example phyloseq object
ps <- readRDS("GG2_physeq.rds")

# Score using default betas
res1 <- compute_risk_score(ps, input_type = "phyloseq")

# Score using custom betas
my_betas <- read.csv("my_betas.csv", stringsAsFactors = FALSE)
res2 <- compute_risk_score(ps, input_type = "phyloseq", betas = my_betas)

# Optional metadata alignment
metadata <- read.table("metadata_agingmicrobiome.txt", header = TRUE, sep = "\t")
res3 <- compute_risk_score(ps, input_type = "phyloseq",
                           metadata = metadata, sampleid_col = "sampleid")
```

---

## Notes

* Make sure to **filter participants before running** if cohort-specific inclusion/exclusion criteria are needed.
* Only genera with sufficient prevalence and abundance contribute to the score.
* CLR transformation is applied automatically.
* The package includes `betas_default` internally for easy use.

---

## License

MIT
