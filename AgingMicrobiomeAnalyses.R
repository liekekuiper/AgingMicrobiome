rm(list=ls())
# Load libraries
library("phyloseq")
library("qiime2R")
library("survival")
library("pcaPP")
library("data.table")
library("vegan")

# Set file paths
modsfile <- "mods_agingmicrobiome.txt"
outsfile <- "outs_agingmicrobiome.txt"
metafile <- "metadata_agingmicrobiome.txt"
subsfile <- "subsets_agingmicrobiome.txt"
featuref <- "feature.table.gg2-2022.10.qza"
taxofile <- "df.gg2.taxonomy.qza"
cohortname <- "cohort.txt"

#Load the functions
source('functions_microbiome_aging.R')

# Load data
models <- scan(modsfile, what = "", sep = "\n")
outcomes <- scan(outsfile, what = "", sep = "\n")
cohort <- scan(cohortname, what = "", sep = "\n")
subsets <- scan(subsfile, what = "", sep = "\n")

#Read in metadata and make sample_data file to be able to add to phyloseq object
metadata <- read.table(metafile, header = TRUE)
rownames(metadata) <- metadata$SampleID
names(metadata)[names(metadata) != "SampleID"] <- tolower(names(metadata)[names(metadata) != "SampleID"])
metadata <- phyloseq::sample_data(metadata)
savenames = tolower(names(metadata))

# Create phyloseq object
physeq_all <- merge_phyloseq(
  qza_to_phyloseq(features = featuref, taxonomy = taxofile),
  metadata
)

print("Phyloseq object created")

# Main loop
for (subs in subsets) {
  physeq <- phy_subset(physeq_all, subs)
  print(paste("Starting analyses of", subs))
  
  physeq <- calculate_alpha_diversity(physeq)
  print("Calculated alphas")
  
  for (method in c("bray", "jaccard")) {
    physeq <- calculate_distance_matrix(physeq, method)
  }
  print("Calculated bray and jaccard")
  physeq_clr <- microbiome::transform(physeq, 'clr')
  print("CLR transformed")
  for (method in c("euclidean", "kendall")) {
  physeq_clr <- calculate_distance_matrix(physeq_clr, method)
  }
  print("Calculated Aitchinson and Kendall")
  
  #Add genera to sample data and assign names
  otu_genus <- as.data.frame(t(otu_table(physeq_clr)))
  colnames(otu_genus) <- as.data.frame(tax_table(physeq_clr))$Genus
  otu_genus <- sample_data(otu_genus)
  physeq_clr <- merge_phyloseq(physeq_clr, otu_genus)
  rm(otu_genus)
  
  #Make datafile ready for analyses.
  datafile <- as.data.frame(sample_data(physeq_clr))
  datafile <- as.data.frame(lapply(datafile, unclass))
  names(datafile) <- tolower(names(datafile))
  
  #Select added variables as determinants
  variables <- setdiff(names(datafile), savenames)
  
  #Prepare dataframe
  results_df <- data.frame(
    Datasplit = character(0),
    Outcome = character(0),
    Variable = character(0),
    Model = character(0),
    N = numeric(0),
    Ncases = numeric(0),
    Coefficient = numeric(0),
    Std.Error = numeric(0),
    HR = numeric(0),
    t.value = numeric(0),
    P = numeric(0)
  )
  print(paste("Performing analyses of", subs))
  for (out in outcomes) {
    for (var in variables) {
      for (model in models) {
        if (!out %in% c("age", "mortality")) {
          model <- paste(model, "age", sep = "+")
        }
        if (grepl("men", subs)) {
          model <- gsub("sex\\+|sex", "", model)
        }
        if (out != "mortality") {
          formula <- as.formula(paste(out, "~", model, "+", var))
          lmmodel <- lm(formula, data = datafile)
          coefficients <- summary(lmmodel)$coefficients
          variable <- var
          result_row <- data.frame(
            Datasplit = subs,
            Outcome = out,
            Variable = var,
            Model = model,
            N = length(lmmodel$fitted.values),
            Ncases = NA,
            Coefficient = coefficients[variable, "Estimate"],
            Std.Error = coefficients[variable, "Std. Error"],
            HR = NA,
            t.value = coefficients[variable, "t value"],
            P = coefficients[variable, "Pr(>|t|)"]
          )
        } else if (out == "mortality") {
          datafile$followup <- datafile$studytime + datafile$age
          formula <- as.formula(paste("Surv(followup,mortality) ~ ", model, "+", var))
          fit <- coxph(formula, data = datafile)
          coefficients <- summary(fit)$coefficients
          variable <- var
          result_row <- data.frame(
            Datasplit = subs,
            Outcome = out,
            Variable = var,
            Model = model,
            N = fit$n,
            Ncases = fit$nevent,
            Coefficient = coefficients[variable, "coef"],
            Std.Error = coefficients[variable, "se(coef)"],
            HR = coefficients[variable, "exp(coef)"],
            t.value = NA,
            P = coefficients[variable, "Pr(>|z|)"]
          )
        }
        results_df <- rbind(results_df, result_row)
      }
    }
  }
  print(paste("Finished analyses of", subs))
  data.table::fwrite(results_df, paste0("Results_", subs, "_", cohort, "_", Sys.Date(), ".csv"))
}
