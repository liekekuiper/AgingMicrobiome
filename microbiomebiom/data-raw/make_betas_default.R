betas_default <- read.csv(
  "data-raw/betas_default.csv",
  stringsAsFactors = FALSE
)

stopifnot(all(c("RISK", "B") %in% colnames(betas_default)))

usethis::use_data(betas_default, overwrite = TRUE)

