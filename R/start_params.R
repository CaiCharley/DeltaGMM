#' Set starting parameters for NLS curve fitting
#' 
#' Construct a set of initial conditions for curve fitting using nonlinear least
#' squares using a mixture of Gaussians. This method finds local maxima within 
#' the chromatogram, orders them by their magnitude from local maxima and uses 
#' the positions and heights of these local maxima (+/- some random noise) as 
#' initial conditions. The "random" method simply picks random values within 
#' the fraction and intensity intervals as starting points for Gaussian 
#' curve-fitting. The initial value of sigma is set by default to a random
#' number within +/- 0.5 of two for both modes; this is based on our manual 
#' inspection of a large number of chromatograms.
#' Derived from PrInCE's `make_initial_conditions`
#'
#' @param chromatogram A numeric vector of a smoothed protein chromatogram
#' @param n_gaussians The number of gaussian parameters to output
#' @param sigma_default The default mean initial value of sigma
#' @inheritParams fit_n_gaussians
#'
#' @importFrom dplyr first nth last
#'
#' @returns A named list containing `A`, `mu`, `sigma` each with length of
#' `n_gaussians`
#'
#' @export
start_params <- function(chromatogram,
                         n_gaussians,
                         gmmctrl = setgmmctrl(),
                         sigma_default = 2) {
  minHeight <- min(chromatogram)
  maxHeight <- max(chromatogram)
  fractions <- length(chromatogram)
  if (gmmctrl$init_method == "guess") {
    peaksX <- which(diff(sign(diff(chromatogram))) == -2) + 1
    if (first(chromatogram) > nth(chromatogram, 2))
      peaksX <- c(1, peaksX)
    if (last(chromatogram) > nth(chromatogram, -2))
      peaksX <- c(peaksX, fractions)
    peaksY <- chromatogram[peaksX]
    peaksX <- peaksX[order(-peaksY)]
    peaksY <- peaksY[order(-peaksY)]

    A <- mu <- sigma <- numeric(0)
    for (i in seq_len(n_gaussians)) {
      A[i] <- peaksY[i] + runif(1) - 0.5
      mu[i] <- peaksX[i] + runif(1, max = 3) - 1.5
      sigma[i] <- sigma_default + runif(1) - 0.5
    }
    A[is.na(A)] <- runif(sum(is.na(A)), min = minHeight, max = maxHeight)
    mu[is.na(mu)] <- runif(sum(is.na(mu)), min = 1, max = fractions)
  } else if (gmmctrl$init_method == "random") {
    A <- runif(n_gaussians, min = minHeight, max = maxHeight)
    mu <- runif(n_gaussians, min = 1, max = fractions)
    sigma <- sigma_default + runif(n_gaussians) - 0.5
  }
  return(list(A = A, mu = mu, sigma = sigma))
}

#' @export
smooth_chromatogram <- function(chromatogram,
                                impute_NA = FALSE,
                                smooth_width = 4) {
  # Benchmark to see if performace matters
  # impute neighbours
  if (impute_NA) {
    fractions <- length(chromatogram)
    for (i in which(is.na(chromatogram))) {
      prevIdx <- i - 1
      nextIdx <- i + 1
      if (prevIdx < 1 | nextIdx > fractions - 1) next
      neighbors <- c(chromatogram[prevIdx], chromatogram[nextIdx])
      if (!any(is.na(neighbors))) {
        chromatogram[i] <- mean(neighbors)
      }
    }
  }

  chromatogram[is.na(chromatogram)] <- 0

  smooth_input <- c(rep(0, smooth_width), chromatogram, rep(0, smooth_width))
  smoothed <- forecast::ma(smooth_input, smooth_width)
  smoothed <- smoothed[(smooth_width + 1):(length(smoothed) - smooth_width)]

  return(smoothed)
}



