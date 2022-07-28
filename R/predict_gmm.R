#' Returns the fitted GMM model from the output of [fit_protein()]
#'
#' @param model a object output from [fit_n_gaussians()] containing `coefs`,
#' `n_gaussians`, and `fractions` from the NLS model
#'
#' @return a vector of the fitted GMM
#' @export
predict_gmm <- function(model) {
  # TODO: cleanup function. Perhaps n_gaussians can be calculated? Fractions not necessary
  coefs <- split(model$coefs, rep(seq_len(3), each = model$n_gaussians))
  coefs <- setNames(coefs, c("A", "mu", "sigma"))

  A <- coefs[["A"]]
  mu <- coefs[["mu"]]
  sigma <- coefs[["sigma"]]

  rowSums(sapply(seq_len(model$n_gaussians),
                 function(i) A[i] * exp(-((seq_len(model$fractions) - mu[i]) / sigma[i])^2)))
}
