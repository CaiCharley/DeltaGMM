% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AICc.R
\name{AICc}
\alias{AICc}
\title{Model selection for Gaussian mixture models}
\usage{
AICc(model)
}
\arguments{
\item{model}{The model fitted by nls to calculate AICc}
}
\value{
the AICc score
}
\description{
Calculate the AIC, corrected AIC, or BIC for a curve fit with a Gaussian
mixture model by nonlinear least squares optimization. This function
permits the calculation of the AIC/AICc/BIC after rejecting some Gaussians
in the model, for example because their centres are outside the bounds of
the profile.
}
