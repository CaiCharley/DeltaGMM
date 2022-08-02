#' Score a fitted model from the output of [fit_protein()] using the F-statistic
#'
#' From the individually fit and simultaneously fit GMMs from chromatograms from
#' two conditions, calculate the F-statistic based score from the differences of
#' RSS values of the fitted models.
#'
#' @param fit A list containing two individually fit models and a null fit model
#' from the output of [fit_protein()]. Models must contain `RSS`, `residual`.
#'
#' @return The protein score
#'
#' @export
score_protein <- function(fit) {
  RSS0 <- fit$null_model$RSS
  RSS1 <- fit$raw1_model$RSS + fit$raw2_model$RSS

  d1 <- fit$null_model$degFreedom
  d2 <- fit$raw1_model$degFreedom + fit$raw2_model$degFreedom

  # d1 <- df.residual(fit$null_model$fit)
  # d2 <- df.residual(fit$raw1_model$fit) + df.residual(fit$raw2_model$fit)

  return((d2 / d1) * (RSS0 / RSS1))
}
