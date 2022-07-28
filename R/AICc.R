#' Model selection for Gaussian mixture models
#'
#' Calculate the AIC, corrected AIC, or BIC for a curve fit with a Gaussian
#' mixture model by nonlinear least squares optimization. This function
#' permits the calculation of the AIC/AICc/BIC after rejecting some Gaussians
#' in the model, for example because their centres are outside the bounds of
#' the profile.
#'
#' @param model The model fitted by nls to calculate AICc
#'
#' @return the AICc score
#'
#' @export

AICc <- function(model) {
  # first, calculate AIC
  AIC <- AIC(model)
  # second, calculate AICc
  # TODO: correct nobs to use passed in parameter as no NAs are passed to nls
  N <- nobs(model) # number of non NA
  k <- length(unlist(coef(model), use.names = FALSE)) + 1
  AICc <- AIC + (2 * k * (k + 1)) / (N - k - 1)
  return(AICc)
}

