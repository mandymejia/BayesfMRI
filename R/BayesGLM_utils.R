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

#' Check design matrix or matrices
#' 
#' @param design,design_type,nS,nT See \code{"BayesGLM"}
#' @return The design matrix as an \code{nS}-length list.
#' @keywords internal 
#' 
BayesGLM_check_design <- function(design, design_type, nS, nT) {
  design_type <- match.arg(
    design_type, c("design", "compare", "per_loc", "onsets")
  )
  if (!is.list(design) || !is.list(design$design)) {
    stop("`design` is not formatted correctly. See `make_design`.")
  }
  if (length(design$design) != nS) {
    stop("The length of `design$design` should match the number of BOLD sessions, ", nS, ".")
  }
  # (Repeat) the checks from `make_design`.
  # `nD`: num. of designs (>1 for per-location or multi modeling). Not used.
  stopifnot(all(vapply(design$design, is.numeric, FALSE)))
  stopifnot(length(unique(lapply(design$design, dim)))==1)
  if (design_type %in% c("regular", "from_onsets")) {
    stopifnot(all(vapply(design$design, is.matrix.or.df, FALSE)))
    # nD <- 1
  } else if (design_type %in% c("multi", "per_location")) {
    stopifnot(all(vapply(design$design, is.array, FALSE)))
    # nD <- dim(design$design[[1]])[3]
  } else { stop() }

  design
}

#' Check nuisance matrix or matrices
#' 
#' @param nuisance,nS,nT See \code{"BayesGLM"}
#' @return The nuisance matrix as an \code{nS}-length list.
#' @keywords internal 
#' 
BayesGLM_check_nuisance <- function(design, design_type, nS, nT) {
  if (is.matrix(nuisance)) { nuisance <- list(single_sess=nuisance) }

  if (!is.null(nuisance)) {
    stopifnot(all(vapply(nuisance, is.matrix, FALSE)))
    stopifnot(all(vapply(nuisance, is.numeric, FALSE)))
  } else {
    nuisance <- vector("list", length = nS)
  }
  nuisance
}

#' Make DCT bases
#' @param hpf,nT,TR See \code{BayesGLM_cifti}
#' @return The matrix of DCT bases, or \code{NULL} if none
#' @importFrom fMRItools dct_bases is_posNum dct_convert
#' @keywords internal 
BayesGLM_cifti_make_DCT <- function(hpf, nT, TR){
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
      nDCT <- round(dct_convert(nT, TR, f=hpf))
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
#' Checks arguments for \code{BayesGLM} and \code{BayesGLM_cifti}
#'
#' Avoids duplicated code between \code{BayesGLM} and \code{BayesGLM_cifti}
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
    Bayes,
    EM,
    ar_order,
    ar_smooth,
    aic,
    n_threads,
    return_INLA,
    verbose,
    meanTol,
    varTol,
    emTol
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
#' @param spde \code{NULL} or the SPDE.
#' @keywords internal
#' @return A list of two: \code{T} for the total number of locations, and 
#'  \code{D} for the number of data locations. If \code{spatial} is provided for
#'  voxel data, there is also \code{DB} for the number of data locations plus 
#'  the number of boundary locations.
get_nV <- function(spatial, type=c("mesh", "voxel"), spde=NULL){
  type <- match.arg(type, c("mesh", "voxel"))

  out <- switch(type, 
    mesh = list(
      T=nrow(spatial$surf$vertices),
      D=sum(spatial$mask)
    ),
    voxel = list(
      T=prod(dim(spatial$label)), # [TO DO] redefine?
      D=sum(spatial$label!=0)
    )
  )

  if (type=="voxel" && !is.null(spde)) {
    out <- c(out, list(DB=spde$n.spde))
  }

  out
}

nT_message <- function(nT) {
  if (length(nT)==1) { return(nT) }
  cat(min(nT), "-", max(nT))
}