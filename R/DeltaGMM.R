#' Scores differential protein interactions in CF-MS data
#'
#' Identifies changes in a protein's complex membership between two conditions
#' by fitting a Gaussian Mixture Model to the protein chromatogram in each
#' condition as a null model where a GMM is fit to both conditions
#' simulatenously. Proteins are then ranked by some metric (such as the
#' F-statistic) to identify differential protein interactions.
#'
#' @param condition1,condition2 A pair of matrices or a paired list of matrices
#'  representing replicates of chromatograms. Proteins are represented as rows
#'  and fractions are represented as columns. Replicates must have matching
#'  number of fractionations. Only proteins appearing in both conditions will be
#'  used.
#'
#' @param score The scoring metric to use. Defaults to `"fstat"`.
#' @param imputectrl An optional named list of class "imputectrl" of parameters
#'  for imputation. See [setimputectrl()] for names of settable control values and
#'  their effects.
#' @param gmmctrl An optional named list of class "gmmctrl" of parameters
#'  for imputation. See [setgmmctrl()] for names of settable control values and
#'  their effects.
#' @param parallel Use multiple cores for fitting.
#'
#' @return A named list of proteins with scores.
#'
#' @export
deltaGMM <- function(condition1,
                     condition2,
                     score = c("fstat", "manhattan"),
                     imputectrl = setimputectrl(),
                     gmmctrl = setgmmctrl(),
                     parallel = TRUE) {
  # parse arguments
  score <- match.arg(score)
  if (!inherits(imputectrl, "imputectrl"))
    stop("controls is not a imputectrl object - see help(\"imputectrl\")")
  if (!inherits(gmmctrl, "gmmctrl"))
    stop("controls is not a gmmctrl object - see help(\"gmmctrl\")")

  # convert conditions to list of matricies
  if (is.matrix(condition1)) {
    condition1 <- list(condition1)
  } else if (!is.list(condition1) || !is.matrix(condition1[[1]])) {
    stop("condition1 must be a matrix or list of matricies")
  }

  if (is.matrix(condition2)) {
    condition2 <- list(condition2)
  } else if (!is.list(condition2) || !is.matrix(condition2[[1]])) {
    stop("condition2 must be a matrix or list of matricies")
  }

  # check all replicates have same number of fractions
  if (length(unique(sapply(c(condition1, condition2), ncol)) != 1)) {
    stop("Replicates must have matching number of fractions")
  }

  # get union of proteins between replicates and intersect between conditions
  condition1_proteins <- purrr::reduce(lapply(condition1, rownames), union)
  condition2_proteins <- purrr::reduce(lapply(condition2, rownames), union)

  proteins <- intersect(condition1_proteins, condition2_proteins)

  # TODO: Proprocessing steps
  # - Set missing values to NA
  # - Use PrInCE::clean_profile() ?

  # TODO: Parallelize fit and scoring
  fit_and_score <- function(protein) {
    raw1 <- DeltaGMM::collect_replicates(protein, condition1)
    raw2 <- DeltaGMM::collect_replicates(protein, condition2)

    # fit each condition individually and both simultanously
    to_fit <- list(raw1_model = raw1,
                   raw2_model = raw2,
                   null_model = rbind(raw1, raw2))

    lapply(to_fit, DeltaGMM::fit_protein,
           imputectrl = imputectrl, gmmctrl = gmmctrl)
  }

  # parallelization
  cl_option <- NULL
  if (parallel) {
    cl_option <- parallel::detectCores()
    if (.Platform$OS.type == "windows") {
      cl_option <- parallel::makeCluster(cl_option)
    }
  }

  protein_fits <- pbapply::pblapply(proteins, fit_and_score, cl = cl_option)

  if (parallel && .Platform$OS.type == "windows") {
    parallel::stopCluster(cl_option)
  }

  return(protein_fits)
}

#' Pull protein's chromatograms across replicates into matrix
#'
#' @param protein The protein name as it appears in the rowname of `replicates`
#' @param replicates A list of matrices representing replicates with protein
#' names as rownames and fractions as columns.
#'
#' @return A matrix with the each replicate of the protein in rows
collect_replicates <- function(protein, replicates) {
  combine <- lapply(replicates, function(rep) {
    rep[rownames(rep) == protein]
  })
  return(purrr::reduce(combine, rbind))
}
