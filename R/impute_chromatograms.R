#' Impute a matrix of chromatograms based on provided parameters
#'
#' @inheritParams fit_n_gaussians
#'
#' @returns a matrix of chromatograms with imputed missing values
#'
#' @export
impute_chromatograms <- function(chromatograms, imputectrl = setimputectrl()) {
  switch(
    imputectrl$imputation,
    "dynamic" = {
      chromatograms[is.na(chromatograms)] <- 0
    },
    "fill_gaps" = {
      chromatograms <- t(apply(chromatograms, 1, function(rep) {
        rep[distance_from_quantified(rep) >= imputectrl$max_dist] <- 0
        return(rep)
      }))
    },
    "weights" = {
      chromatograms[is.na(chromatograms)] <- 0
    },
    "none" = { },
  )

  return(chromatograms)
}

#' Find distance to nearest quantified fraction
#'
#' @param chrom A numeric vector of a protein chromatogram
#'
#' @returns A numeric vector with same length as chrom with distance from
#' nearest quantified fraction for each fraction
distance_from_quantified <- function(chrom) {
  return(
    pmin(
      distance_from_left(chrom),
      rev(distance_from_left(rev(chrom)))
    )
  )
}

distance_from_left <- function(chrom) {
  distances <- numeric(length = length(chrom))
  distance <- length(chrom)
  for (fraction in seq_along(chrom)) {
    if (is.null(chrom[fraction]) |
      is.na(chrom[fraction]) |
      is.nan(chrom[fraction])) {
      distance <- distance + 1
    } else {
      distance <- 0
    }
    distances[fraction] <- distance
  }
  return(distances)
}

