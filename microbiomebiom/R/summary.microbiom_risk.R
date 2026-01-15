#' Summary method for microbiomebiom risk scores
#'
#' @param object An object of class microbiomebiom_risk
#' @param ... Not used
#' @export
summary.microbiomebiom_risk <- function(object, ...) {

  cat("Microbiomebiom risk score summary\n")
  cat("---------------------------------\n")
  cat("Number of samples:      ", object$n_samples, "\n")
  cat("Number of genera used:  ", object$n_genera_used, "\n")

  invisible(object)
}

