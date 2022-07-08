#' Fits a Gaussian Mixture model to protein chromatograms
#'
#' Fits one GMM to the raw chromatograms in each condition and one GMM too all
#' the chromatograms as if they were from the same condition to make a null
#' model. A specified score is then calculated from features of the fits to rank
#' the amount of differential change
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
                        imputectrl = imputectrl(),
                        gmmctrl = gmmctrl()) {
  if (!is.matrix(raws) && !is.vector(raws)) {
    stop("raws must be matrix or single vector")
  } else if (is.vector(raws)) {
    raws <- matrix(raws, nrow = 1)
  }

  # find number of fractions that are quantified in at least one replicate
  quantified_fractions <- length(unique(which(!is.na(raws), arr.ind = T)[, 2]))

  # don't fit mixtures with more parameters than (experimental) points
  if (!is.null(quantified_fractions)) {
    max_gaussians <- min(max_gaussians, floor(quantified_fractions / 3))
  }

  # fit models
  fits <- list()
  for (n in seq_len(gmmctrl$max_gaussians)) {
    fits[[n]] <- DeltaGMM::fit_n_gaussians(chormatograms,
                                           n_gaussian = n,
                                           imputectrl = imputectrl,
                                           gmmctrl = gmmctrl)
  }

  # choose model
  if (length(fits) == 0) {
    return(NULL)
  } else {
    # remove any models that failed to fit
    drop <- purrr::map_lgl(fits, is.null)
    fits <- fits[!drop]

    criterions <- lapply(fits, function(x) criterion_fn(x$fit))
    # return NULL if all criteria are NA
    if (sum(is.na(criterions)) == length(criterions)) return(NULL)

    return(fits[[which.min(criterions)]])
  }
}
