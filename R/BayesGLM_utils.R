#' Scale the design matrix
#'
#' @param design_mat The original (unscaled) design matrix that is T x K, where
#'     T is the number of time points, and k is the number of field covariates
#'
#' @return A scaled design matrix
#'
#' @keywords internal
#'
scale_design_mat <- function(design_mat) {
  stopifnot(is.matrix(design_mat))
  apply(design_mat,2,function(field) {
    #(field - mean(field)) / max(abs(field))
    field / max((field))
  })
}

#' Make DCT bases
#' @param hpf,nT,TR See \code{BayesGLM}
#' @return The matrix of DCT bases, or \code{NULL} if none
#' @importFrom fMRItools dct_bases is_posNum dct_convert
#' @keywords internal
BayesGLM_make_DCT <- function(hpf, nT, TR){
  DCT <- NULL # used to be an argument.

  if (!is.null(DCT)) {
    stopifnot(is_posNum(DCT, zero_ok=TRUE) && DCT==round(DCT))
    if (DCT==0) { DCT <- NULL }
  }
  if (!is.null(hpf)) {
    stopifnot(is_posNum(hpf, zero_ok=TRUE))
    if (hpf==0) { hpf <- NULL }
  }

  if (!is.null(hpf) || !is.null(DCT)) {
    # Get the num. of bases for this session.
    if (!is.null(hpf)) {
      nDCT <- round(dct_convert(T_=nT, TR=TR, f=hpf))
    } else {
      nDCT <- DCT
    }
    if (nDCT > 0) {
      fMRItools::dct_bases(nT, nDCT)
    } else {
      NULL
    }
  } else {
    NULL
  }
}

#' Bayes GLM arg checks
#'
#' Checks arguments for \code{BayesGLM} and \code{fit_bayesglm}
#'
#' Avoids duplicated code between \code{BayesGLM} and \code{fit_bayesglm}
#'
#' @param scale_BOLD See \code{\link{BayesGLM}}.
#' @param Bayes,EM See \code{\link{BayesGLM}}.
#' @param ar_order,ar_smooth,aic See \code{\link{BayesGLM}}.
#' @param n_threads See \code{\link{BayesGLM}}.
#' @param return_INLA See \code{\link{BayesGLM}}.
#' @param verbose See \code{\link{BayesGLM}}.
# @param combine_sessions See \code{\link{BayesGLM}}.
#' @param meanTol,varTol,emTol See \code{\link{BayesGLM}}.
#'
#' @return The arguments that may have changed, in a list: \code{scale_BOLD},
#'  \code{do_Bayesian}, \code{do_EM}, and \code{do_pw}.
#'
#' @importFrom fMRItools is_1
#' @keywords internal
BayesGLM_argChecks <- function(
    scale_BOLD=c("mean", "sd", "none"),
    Bayes=TRUE,
    EM=FALSE,
    ar_order=6,
    ar_smooth=5,
    aic=FALSE,
    n_threads=4,
    return_INLA=c("trimmed", "full", "minimal"),
    verbose=1,
    meanTol=1e-6,
    varTol=1e-6,
    emTol=1e-3
){

  if (isTRUE(scale_BOLD)) {
    message("Setting `scale_BOLD` to 'mean'"); scale_BOLD <- "mean"
  }
  if (isFALSE(scale_BOLD)) {
    message("Setting `scale_BOLD` to 'none'"); scale_BOLD <- "none"
  }
  scale_BOLD <- match.arg(scale_BOLD, c("mean", "sd", "none"))

  stopifnot(fMRItools::is_1(Bayes, "logical"))
  stopifnot(fMRItools::is_1(EM, "logical"))
  if (EM && !Bayes) {
    warning("EM is a Bayesian method: setting `Bayes` to `TRUE`.")
    Bayes <- TRUE
  }
  if (Bayes) {
    if (!EM) { check_INLA(require_PARDISO=FALSE) }
  }

  if(EM) stop("EM not available.") #not currently available

  if (is.null(ar_order)) ar_order <- 0
  stopifnot(fMRItools::is_1(ar_order, "numeric"))
  do_pw <- ar_order > 0
  if (is.null(ar_smooth)) ar_smooth <- 0
  stopifnot(fMRItools::is_1(ar_smooth, "numeric"))
  stopifnot(fMRItools::is_1(aic, "logical"))

  stopifnot(fMRItools::is_1(n_threads, "numeric"))
  stopifnot(n_threads <= parallel::detectCores())

  if (isTRUE(return_INLA)) {
    message("Setting `return_INLA` to 'trimmed'"); return_INLA <- "trimmed"
  }
  if (isFALSE(return_INLA)) {
    message("Setting `return_INLA` to 'minimal'"); return_INLA <- "minimal"
  }
  return_INLA <- match.arg(return_INLA, c("trimmed", "full", "minimal"))

  if (isTRUE(verbose)) { verbose <- 2 }
  if (isFALSE(verbose)) { verbose <- 0 }
  stopifnot(fMRItools::is_posNum(verbose, zero_ok=TRUE))

  stopifnot(fMRItools::is_posNum(meanTol))
  stopifnot(fMRItools::is_posNum(varTol))
  stopifnot(fMRItools::is_posNum(emTol))

  # Return new parameters, and parameters that may have changed.
  list(
    scale_BOLD=scale_BOLD,
    Bayes=Bayes,
    EM = EM,
    do_pw = do_pw,
    return_INLA=return_INLA
  )
}

#' Get number of locations for various masks
#'
#' Get number of locations for various masks: total, model, data.
#'
#' @param spatial \code{spatial}
#' @keywords internal
#' @return A list of two: \code{T} for the total number of locations, and
#'  \code{D} for the number of data locations. If \code{spatial} is provided for
#'  voxel data, there is also \code{DB} for the number of data locations plus
#'  the number of boundary locations.
get_nV <- function(spatial){

  # for `nV_total`, use these instead of masks, b/c `BayesGLM`
  #   will have these before the masks are read in.
  nV_total <- switch(spatial$spatial_type,
    vertex = nrow(spatial$surf$vertices),
    voxel = prod(dim(spatial$maskIn))
  )

  nV_input <- sum(spatial$maskIn)
  nV_model <- spatial$nV_M
  nV_mdata <- sum(spatial$maskMdat)
  nV_mbuffer <- nV_model - sum(spatial$maskMdat)

  list(
    total = nV_total, #length(spatial$maskIn),
    input = nV_input,
    model = nV_model,
    mdata = nV_mdata,
    mbuffer = nV_mbuffer
  )
}

nT_message <- function(nT) {
  if (length(nT)==1) { return(nT) }
  cat(min(nT), "-", max(nT))
}

#' Set column values to zero for sparse matrix
#'
#' Set column values to zero for sparse matrix
#' @param mat dgCMatrix
#' @param cols The column indices to set to zero
#' @return The modified sparse matrix
#' @keywords internal
dgCMatrix_cols_to_zero <- function(mat, cols) {
  stopifnot(inherits(mat, "dgCMatrix"))
  stopifnot(is.numeric(cols))

  # column of each nonempty value in the sparse matrix: https://stackoverflow.com/questions/21099612/extract-i-and-j-from-a-sparse-matrix
  cols_existing <- findInterval(seq(mat@x)-1,mat@p[-1])+1

  # update i and x
  # meaning of i, x, p: https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
  # delete values in the columns to drop: https://stackoverflow.com/questions/33775291/r-matrix-set-particular-elements-of-sparse-matrix-to-zero
  sparseVals_mask <- !(cols_existing %in% cols)
  new_i <- mat@i[sparseVals_mask]
  new_x <- mat@x[sparseVals_mask]

  # now to update p:
  # https://stackoverflow.com/questions/20008200/r-constructing-sparse-matrix
  # "the difference between the ith and the (i-1)th element in p is
  # the number of x elements in (i-1) column"
  # new_p <- mat@p
  nv_cols <- diff(mat@p) # num elements in each column
  nv_cols[cols] <- 0 # set to zero where needed
  new_p <- c(0, cumsum(nv_cols)) # go back to the cumulative sum form

  mat@x <- new_x
  mat@i <- new_i
  mat@p <- as.integer(new_p)

  mat
}

#' Validate \code{spatial}
#'
#' Validate \code{spatial}
#'
#' @param spatial \code{spatial}
#' @return \code{NULL}, invisibly
#' @keywords internal
#'
validate_spatial <- function(spatial) {
  stopifnot(is.list(spatial))
  if (spatial$spatial_type == "vertex") {
    stopifnot(length(spatial) == 7)
    stopifnot(names(spatial) == c(
      "spatial_type", "surf", "maskIn", "maskMdat", "nV_M", "Mmap", "mesh"
    ))
  } else if (spatial$spatial_type == "voxel") {
    stopifnot(length(spatial) == 11)
    stopifnot(names(spatial) == c(
      "spatial_type", "labels", "trans_mat", "trans_units",
      "nbhd_order", "buffer",
      "maskIn", "labsMdat", "maskMdat", "nV_M", "Mmap")
    )
  } else { stop("Unknown spatial$spatial_type.") }

  invisible(NULL)
}
