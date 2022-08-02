#' Fits a Gaussian Mixture model to protein chromatograms
#'
#' Attempts to fit multiple GMM to the provided chromatograms with different
#' number of gaussian components from 1 to `max_gaussians`. Selects best model
#' based on the specified information criterion in `gmmctrl`.
#'
#' @param chromatograms A matrix of raw chromatograms in rows (replicates)
#' and fractions in columns to fit gaussians to. A single numeric vector will be
#' coerced to a matrix with one row.
#' @inheritParams deltaGMM
#'
#' @return A named list containing the fitted model for the protein and other
#' details
#'
#' @export
fit_protein <- function(chromatograms,
                        imputectrl = setimputectrl(),
                        gmmctrl = setgmmctrl()) {
  # parse args
  if (!is.matrix(chromatograms) && !is.vector(chromatograms)) {
    stop("chromatograms must be matrix or single vector")
  } else if (is.vector(chromatograms)) {
    chromatograms <- matrix(chromatograms, nrow = 1)
  }

  # find number of fractions that are quantified in at least one replicate
  quantified_fractions <- sum(!(colSums(!is.na(chromatograms)) == 0))

  # don't fit mixtures with more parameters than (experimental) points
  if (!is.null(quantified_fractions)) {
    gmmctrl$max_gaussians <-
      min(gmmctrl$max_gaussians, floor(quantified_fractions / 3))
  }

  # fit models
  fits <- lapply(seq_len(gmmctrl$max_gaussians), function(n) {
    fit_n_gaussians(chromatograms,
                         n_gaussians = n,
                         imputectrl = imputectrl,
                         gmmctrl = gmmctrl)
  })

  # choose model
  if (length(fits) == 0) {
    return(NULL)
  } else {
    # remove any models that failed to fit
    fits <- fits[!sapply(fits, is.null)]

    criterions <- sapply(fits, function(x) x$criterion)

    # return NULL if all criteria are NULL
    if (all(sapply(criterions, is.null)))
      return(NULL)
    else
      return(fits[[which.min(criterions)]])
  }
}
