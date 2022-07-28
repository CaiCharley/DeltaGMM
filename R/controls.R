#' Parameters used to control imputation of chromatograms
#'
#' @param imputation Type of imputation to use for missing values.
#'  * `"dynamic"` (the default): Imputes NAs as 0, then dynamically sets the
#'  weight based on the distance to the nearest quantified fraction with
#'  [distance_from_quantified()]. Weights are passed to [nls()].
#'
#'  * `"weights"`: Used with `na_weights` argument. Sets all missing values to 0
#'  and sets a fixed weight based on `na_weights`.
#'
#'  * `"fill_gaps"`: Used with `max_dist` argument. Sets missing values that are
#'  greater or equal to`max_dist` from the nearest quantified fraction to 0 and
#'  leaves closer missing values as NA. This was the first method of imputation
#'  developed to make sure nls doesn't fit gaussians to large gaps without
#'  penalization.
#'
#'  * `"none"`: No imputation and raw values are passed to [nls()].
#' @param max_dist Maximum gap of unquantified fractions before imputation
#' @param na_weight The weight of missing values imputed as 0
#' @param log_k,log_x The `k` and `x_0` value of the logistic function used in
#'  dynamic imputation of weights. `k` represents the value of the sigmoid
#'  midpoint or how many fractions away from the nearest quantified fraction
#'  for the imputed weights to be 0.5. `k` is the steepness of curve transition
#'  from 0 to 1. Default of 1.5 and 4 were determined emperically.
#'
#' @return An named list of class `imputectrl` with specified parameters
#'
#' @export
setimputectrl <- function(imputation = c("dynamic", "fill_gaps",
                                         "weights", "none"),
                          max_dist = 3,
                          na_weight = 0.1,
                          log_k = 1.5,
                          log_x = 4) {
  imputation <- match.arg(imputation)

  if (max_dist < 0)
    stop("max_dist must be positive")
  if (na_weight < 0 || na_weight > 1)
    stop("na_weight must be in range [0, 1]")
  if (log_k < 0 || log_x < 0)
    stop("log_k and log_x must be positive")

  control <- list(imputation = imputation,
                  max_dist = max_dist,
                  na_weight = na_weight,
                  log_k = log_k,
                  log_x = log_x
  )
  class(control) <- "imputectrl"
  return(control)
}

#' Parameters used to control GMM fitting
#'
#' @param max_gaussians The maximum number of Gaussians to fit. Defaults to 5.
#' @param criterion_fn function used in model selection.
#'  Choose one of [AICc()] (default), [AIC()], [BIC()].
#' @param max_iterations The number of times to try fitting the curve with
#' different initial conditions.
#' @param init_method The method used to select the initial conditions for
#'  nonlinear least squares optimization. One of `"guess"` (default) or `
#'  "random"`.
#' @param height_min Ratio to the max value within the chromatograms for which
#'  gaussians with lower heights will be filtered out. Defaults to 0.15.
#' @param height_max Ratio to the max value within the chromatograms for
#'  which gaussians with larger heights will be filtered out. Defaults to 1.5.
#' @param variance_min Minimum variance for fitted gaussians. Defaults to 1.
#' @param variance_max Maximum variance for fitted gaussians. `0` sets maximum
#'  variance to half the fractions (default).
#' @param rssweights Calculate RSS with imputed values adjusted by their weight.
#' @param return_fit Returns the fitted model from [nls()]. Useful for
#' exporatory analysis. Defaults to `FALSE`.
#'
#' @return An named list of class `gmmctrl` with specified parameters
#'
#' @export
setgmmctrl <- function(max_gaussians = 5,
                       criterion_fn = AICc,
                       max_iterations = 5,
                       init_method = c("guess", "random"),
                       height_min = 0.1,
                       height_max = 1.1,
                       variance_min = 2,
                       variance_max = 0,
                       rssweights = TRUE,
                       return_fit = FALSE) {
  init_method <- match.arg(init_method)

  if (max_gaussians < 1)
    stop("max_gaussians must be greater than 0")
  if (max_iterations < 1)
    stop("max_iterations must be greater than 0")
  if (height_min < 0 ||
    height_max < 0 ||
    height_min >= height_max)
    stop("min and max height must be greater than 0 and max > min")
  if (variance_min < 0 ||
    variance_max < 0 ||
    (variance_max != 0 && variance_min >= variance_max))
    stop("min and max variance must be greater than 0 and max > min")

  control <- list(max_gaussians = max_gaussians,
                  criterion_fn = criterion_fn,
                  max_iterations = max_iterations,
                  init_method = init_method,
                  height_min = height_min,
                  height_max = height_max,
                  variance_min = variance_min,
                  variance_max = variance_max,
                  rssweights = rssweights,
                  return_fit = return_fit)
  class(control) <- "gmmctrl"
  return(control)
}
