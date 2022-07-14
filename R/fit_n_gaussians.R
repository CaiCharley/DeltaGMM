#' Fits a Gaussian Mixture Model to chromatograms
#'
#' Attempts to fit a Gaussian Mixture Model to the provided chromatograms with
#' `n_gaussians` number of components using [nls()]. Performs imputation of
#' chromatograms based on `imputectrl` parameters.
#'
#' @inheritParams fit_protein
#' @param n_gaussians The number of gaussian components to fit.
#'
#' @import stats
#'
#' @returns A fitted model of the chromatograms with various details
#'
#' @export
fit_n_gaussians <- function(chromatograms,
                            n_gaussians,
                            imputectrl = setimputectrl(),
                            gmmctrl = setgmmctrl()) {
  # parse args
  if (!is.matrix(chromatograms) && !is.vector(chromatograms)) {
    stop("chromatograms must be matrix or single vector")
  } else if (is.vector(chromatograms)) {
    chromatograms <- matrix(chromatograms, nrow = 1)
  }

  fractions <- ncol(chromatograms)
  replicates <- nrow(chromatograms)

  # set max gaussian variance based on number of fractions
  if (gmmctrl$variance_max == 0) {
    gmmctrl$variance_max <- fractions / 2
  }

  # imputation
  imputed <- impute_chromatograms(chromatograms, imputectrl)
  weights <- impute_weights(chromatograms, imputectrl)

  # join both conditions and with fraction label
  tidy_chromatograms <- data.frame(
    fraction = rep(seq_len(fractions), times = replicates),
    value = as.vector(t(imputed))
  )

  # create smoothed average of both conditions to build initial conditions
  smoothed_profile <- smooth_chromatogram(colMeans(chromatograms, na.rm = T))

  # set max and min values for model.
  # (A) is proportional to highest fraction.
  # (mu) is based on first and last quantified fraction in both conditions.
  # (sigma) is based on parameters
  min_values <- c(A = max(chromatograms, na.rm = T) * gmmctrl$height_min,
                  mu = min(which(!is.na(chromatograms), arr.ind = T)[, 2]),
                  sigma = gmmctrl$variance_min)
  max_values <- c(A = max(chromatograms, na.rm = T) * gmmctrl$height_max,
                  mu = max(which(!is.na(chromatograms), arr.ind = T)[, 2]),
                  sigma = gmmctrl$variance_max)

  bestRSS <- Inf
  bestFit <- NULL
  for (i in seq(gmmctrl$max_iterations)) {

    # make initial conditions with condition average
    # TODO: Move making starting params out of loop, just add noise in loop or
    # using "random"
    initial_conditions <-
      start_params(smoothed_profile, n_gaussians, gmmctrl)

    # trim initial conditions within bounds of min and max values
    # TODO: Move into start_params()?
    initial_conditions <-
      purrr::pmap(list(initial_conditions, min_values, max_values),
                  function(parameter_list, minimum, maximum) {
                    parameter_list[parameter_list < minimum] <- minimum
                    parameter_list[parameter_list > maximum] <- maximum
                    return(parameter_list)
                  })

    A <- initial_conditions$A
    mu <- initial_conditions$mu
    sigma <- initial_conditions$sigma

    # fit the model
    fit <- tryCatch({
      suppressWarnings(
        nls(value ~ gmm_model(fraction, n_gaussians, A, mu, sigma),
            data = tidy_chromatograms,
            start = list(A = A, mu = mu, sigma = sigma),
            algorithm = "port",
            trace = FALSE,
            weights = weights,
            na.action = na.exclude,
            lower = rep(min_values, each = n_gaussians),
            upper = rep(max_values, each = n_gaussians),
            control = list(warnOnly = TRUE, minFactor = 1 / 2048)))
    }, error = function(e) {
      e
    }, simpleError = function(e) {
      e
    })
    if ("error" %in% class(fit))
      next

    # TODO: account for imputed NAs in RSS, resid doesn't account weights, will this affect aic?
    if (gmmctrl$rssweights == T) {
      RSS <- sum((resid(fit) * weights)^2, na.rm = T)
    } else {
      RSS <- sum((resid(fit))^2, na.rm = T)
    }
    # replace best fit with this model?
    if (RSS < bestRSS) {
      bestRSS <- RSS
      bestFit <- fit
    }
  }

  # return bestFit if fitted, otherwise return null
  if (is.null(bestFit)) return(NULL)

  # TODO: add criterion, coefs, df.residual, fix hardcoded AIC
  res <- list(RSS = bestRSS,
              weights = weights,
              n_gaussians = n_gaussians,
              chromatograms = nrow(chromatograms),
              criterion = AIC(bestFit))

  if (gmmctrl$return_fit)
    return(c(bestFit, res))
  else
    return(res)
}

# gaussian mixture model
gmm_model <- function(x, n, m, A, mu, sigma) {
  .colSums(
    A * exp(-(((matrix(rep(x, each = n), nrow = n) - mu) / sigma)^2)),
    n, m
  )
}