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