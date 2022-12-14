% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_protein.R
\name{score_protein}
\alias{score_protein}
\title{Score a fitted model from the output of \code{\link[=fit_protein]{fit_protein()}} using the F-statistic}
\usage{
score_protein(fit)
}
\arguments{
\item{fit}{A list containing two individually fit models and a null fit model
from the output of \code{\link[=fit_protein]{fit_protein()}}. Models must contain \code{RSS}, \code{residual}.}
}
\value{
The protein score
}
\description{
From the individually fit and simultaneously fit GMMs from chromatograms from
two conditions, calculate the F-statistic based score from the differences of
RSS values of the fitted models.
}
