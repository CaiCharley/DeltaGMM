#' Calculate weights for each value in the matrix for [nls()] fitting
#'
#' @inheritParams fit_n_gaussians
#'
#' @returns a matrix of chromatograms with imputed missing values
#'
#' @export
impute_weights <- function(chromatograms, imputectrl = setimputectrl()) {
  # Caution to transpose order of weights for replicates correctly
  switch(
    imputectrl$imputation,
    "dynamic" = {
      distances <- apply(chromatograms, 1, distance_from_quantified)
      weights <-
        1 / (1 + exp(-imputectrl$log_k * (distances - imputectrl$log_x)))
      weights[distances == 0] <- 1
      weights <- as.numeric(weights)
    },
    "weights" = {
      na_indexes <- as.logical(t(is.na(chromatograms)))
      weights <- rep(1, length(chromatograms))
      weights[na_indexes] <- imputectrl$na_weight
    },
    {
      weights <- rep(1, length(chromatograms))
    }
  )

  return(weights)
}
