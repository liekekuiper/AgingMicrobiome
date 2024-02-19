# Create load function
SmartLoad=function(file,newname){
  loadname=load(file)
  return(eval(parse(text=loadname)))
}

# Subset function
phy_subset <- function(physeq, subinfo) {
  if (subinfo == "all") {
    df <- physeq
  } else {
    df <- switch(
      subinfo,
      men = subset_samples(physeq, sex == "men"),
      women = subset_samples(physeq, sex == "women"),
      age_1 = subset_samples(physeq, age >= 18 & age < 40),
      age_2 = subset_samples(physeq, age >= 40 & age < 50),
      age_3 = subset_samples(physeq, age >= 50 & age < 60),
      age_4 = subset_samples(physeq, age >= 60 & age < 70),
      age_5 = subset_samples(physeq, age >= 70)
    )
  }
  
  df <- tax_glom(df, "Genus")
  
  factor_names <- if (subinfo %in% c("women", "men")) {
    c("ppump", "metfor", "statin", "race")
  } else {
    c("sex", "ppump", "metfor", "statin", "race")
  }
  
  sample_data(df)[, factor_names] <- lapply(sample_data(df)[, factor_names], factor)
  return(df)
}

# Alpha diversity function
calculate_alpha_diversity <- function(physeq) {
  alpha_div <- estimate_richness(physeq, split = TRUE, measures = c("Chao1", "Shannon", "InvSimpson"))
  sample_data(physeq)$Chao1 <- alpha_div$Chao1
  sample_data(physeq)$Shannon <- alpha_div$Shannon
  sample_data(physeq)$InvSimpson <- alpha_div$InvSimpson
  return(physeq)
}

# Distance matrix function
calculate_distance_matrix <- function(physeq, method) {
  if(!method %in% c("euclidean", "kendall")){
    distance_matrix <- phyloseq::distance(physeq, method = method, type = "samples")
    df <- as.data.frame(as.matrix(distance_matrix))
    df[df==0]<-NA
    }
  else if(method == "euclidean"){
    distance_matrix <- as.matrix(vegdist(t(otu_table(physeq)),method=method))
    diag(distance_matrix) = NA
    df<-as.data.frame(distance_matrix)
  }else if(method == "kendall"){
    distance_matrix <- as.matrix(1-cor.fk(otu_table(physeq))/2)
    diag(distance_matrix) = NA
    df<-as.data.frame(distance_matrix)
  }
  sample_data(physeq)[,paste0("min_", method)] <- sapply(df, min, na.rm = TRUE)
  return(physeq)
}
