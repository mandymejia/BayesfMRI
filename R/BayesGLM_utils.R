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
#' Checks arguments for \code{BayesGLM} and \code{BayesGLM_fun}
#'
#' Avoids duplicated code between \code{BayesGLM} and \code{BayesGLM_fun}
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
    scale_BOLD,
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
    message("Setting `scale_BOLD` to 'auto'"); scale_BOLD <- "auto"
  }
  if (isFALSE(scale_BOLD)) {
    message("Setting `scale_BOLD` to 'none'"); scale_BOLD <- "none"
  }
  scale_BOLD <- match.arg(scale_BOLD, c("auto", "mean", "sd", "none"))

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
#' @param type \code{"mesh"}, or \code{"voxel"}.
#' @keywords internal
#' @return A list of two: \code{T} for the total number of locations, and
#'  \code{D} for the number of data locations. If \code{spatial} is provided for
#'  voxel data, there is also \code{DB} for the number of data locations plus
#'  the number of boundary locations.
get_nV <- function(spatial, type=c("mesh", "voxel")){
  type <- match.arg(type, c("mesh", "voxel"))

  out <- switch(type,
    mesh = list(
      T=nrow(spatial$surf$vertices),
      D=sum(spatial$mask)
    ),
    voxel = list(
      T=prod(dim(spatial$labels)), # [TO DO] redefine?
      D=sum(spatial$labels!=0)
    )
  )

  if (type=="voxel" && !is.null(spatial$buffer_mask)) {
    stopifnot(out$D == sum(spatial$buffer_mask))
    out <- c(out, list(DB=length(spatial$buffer_mask)))
  }

  out
}

nT_message <- function(nT) {
  if (length(nT)==1) { return(nT) }
  cat(min(nT), "-", max(nT))
}
