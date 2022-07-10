#' Fits a Gaussian Mixture Model to chromatograms
#'
#' Attempts to fit a Gaussian Mixture Model to the provided chromatograms with
#' `n_gaussians` number of components using [nls()]. Performs imputation of
#' chromatograms based on `imputectrl` parameters.
#'
#' @inheritParams fit_protein
#' @param n_gaussians The number of gaussian components to fit.
#'
#' @returns A fitted model of the chromatograms with various details
#'
#' @export
fit_n_gaussians <- function(chromatograms,
                            n_gaussians,
                            imputectrl = imputectrl(),
                            gmmctrl = gmmctrl()) {
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
  # TODO: Can we put in imputed values and replace with simpler smoother?
  smoothed_profile <- PrInCE::clean_profile(
    colMeans(chromatograms, na.rm = T), noise_floor = 0
  )

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
  for (i in seq(gmmctrl$iterations)) {
    # TODO: Replace with fixed version
    # make initial conditions with condition average
    initial_conditions <-
      start_params(smoothed_profile, n_gaussians, gmmctrl$method)

    # trim initial conditions within bounds of min and max values
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
    if (rssweights == T) {
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

  # TODO: add criterion, coefs, df.residual
  res <- list(RSS = bestRSS,
              weights = weights,
              n_gaussians = n_gaussians,
              chromatograms = nrow(raws))

  if (gmmctrl$return_fit)
    return(c(bestFit, res))
  else
    return(res)
}

# gaussian mixture model
gmm_model <- function(x, n, A, mu, sigma) {
  rowSums(sapply(seq_len(n),
                 function(i) A[i] * exp(-((x - mu[i]) / sigma[i])^2)))
}