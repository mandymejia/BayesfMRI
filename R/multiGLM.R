#' multiGLM
#'
#' Performs classical GLM for task fMRI activation, comparing multiple designs
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param BOLD,design,nuisance Session-length list of numeric matrices/arrays,
#'  each with volumes along the first dimension.
#' @param spatial Gives the spatial information. For surface data, a list of two:
#'  \describe{
#'    \item{surf}{A list of two: \code{vertices} \eqn{V \times 3} numeric matrix of vertex locations in XYZ coordinate space, and \code{faces}, \eqn{F \times 3} matrix of positive integers defining the triangular faces.}
#'    \item{mask}{Mask of locations with valid data.}
#' }
#'  For voxel data, a list of six:
#'  \describe{
#'    \item{label}{3D array of labeled locations to include in the model.}
#'    \item{trans_mat}{Projection matrix to convert voxel indices to XYZ position. Can be \code{NULL}.}
#'    \item{trans_units}{XYZ units. Can be \code{NULL}.}
#'    \item{nbhd_order}{See documentation for \code{\link{BayesGLM_cifti}}.}
#'    \item{buffer}{See documentation for \code{\link{BayesGLM_cifti}}.}
#' }
#' @inheritParams session_names_Param
#' @inheritParams scale_BOLD_Param
#' @inheritParams Bayes_Param
# @inheritParams EM_Param
#' @inheritParams ar_order_Param
#' @inheritParams ar_smooth_Param
#' @inheritParams aic_Param
#' @inheritParams n_threads_Param
#' @inheritParams return_INLA_Param
#' @inheritParams verbose_Param
# @inheritParams combine_sessions_Param
#' @param meanTol,varTol Tolerance for mean, variance and SNR of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Default: \code{1e-6} for mean and variance, \code{50} for SNR.
# Note: \code{snrTol} currently not in use, but SNR maps are returned for visualization.
# @inheritParams emTol_Param
#'
#' @return A \code{"CompareGLM"} object: a list with elements
#'  \describe{
#'    \item{field_estimates}{The estimated coefficients for the Bayesian model.}
#'    \item{mask}{A mask of \code{mesh_orig} indicating the locations inside \code{mesh}.}
#'    \item{design}{The design matrix, after centering and scaling, but before any nuisance regression or prewhitening.}
#'    \item{field_names}{The names of the fields.}
#'    \item{session_names}{The names of the sessions.}
#'    \item{hyperpar_posteriors}{Hyperparameter posterior densities.}
#'    \item{theta_estimates}{Theta estimates from the Bayesian model.}
#'    \item{posterior_Sig_inv}{For joint group modeling.}
#'    \item{mu_theta}{For joint group modeling.}
#'    \item{Q_theta}{For joint group modeling.}
#'    \item{y}{For joint group modeling: The BOLD data after any centering, scaling, nuisance regression, or prewhitening.}
#'    \item{X}{For joint group modeling: The design matrix after any centering, scaling, nuisance regression, or prewhitening.}
#'    \item{prewhiten_info}{Vectors of values across locations: \code{phi} (AR coefficients averaged across sessions), \code{sigma_sq} (residual variance averaged across sessions), and AIC (the maximum across sessions).}
#'    \item{call}{match.call() for this function call.}
#'  }
#'
#' @importFrom matrixStats colVars
#' @importFrom Matrix bandSparse bdiag crossprod solve Diagonal
#' @importFrom parallel detectCores makeCluster clusterMap stopCluster
#' @importFrom stats as.formula var
#' @importFrom fMRItools is_1 nuisance_regression scale_timeseries
#'
#' @importFrom utils tail
#'
#' @importFrom methods as
#' @export
multiGLM <- function(
  BOLD,
  design,
  nuisance=NULL,
  # Below arguments shared with `BayesGLM_cifti`.
  scale_BOLD = c("auto", "mean", "sd", "none"),
  verbose = 1,
  meanTol = 1e-6,
  varTol = 1e-6#,
  #snrTol = 50,
  #emTol = 1e-3
  ){

  EM <- FALSE
  emTol <- 1e-3

  # Argument checks. -----------------------------------------------------------
  ### Simple parameters. -------------------------------------------------------
  do <- vector("list")

  # In a separate function because these checks are shared with `BayesGLM`.
  x <- BayesGLM_argChecks(
    scale_BOLD = scale_BOLD,
    Bayes = Bayes,
    EM = EM,
    ar_order = ar_order,
    ar_smooth = ar_smooth,
    aic = aic,
    n_threads = n_threads,
    return_INLA = return_INLA,
    verbose = verbose,
    meanTol = meanTol,
    varTol = varTol,
    emTol = emTol
  )
  scale_BOLD <- x$scale_BOLD
  rm(x)
  rm(EM, emTol) # not used here.

  # Modeled after `BayesGLM_cifti` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Check `BOLD`. ------------------------------------------------------------
  nS <- length(BOLD)
  if(nS!=1) stop("Not supported: multi-session in `compareGLM`.")
  ### Check `design`. ----------------------------------------------------------
  # Make `design` a sessions-length list of design matrices.
  #   Get `nK`, `field_names`, and `do$perLocDesign`. Check for consistent dims
  #   across sessions.
  x <- BayesGLM_format_design(design, nS_expect=nS)
  design <- x$design
  nT <- x$nT
  nK <- x$nK
  nD <- x$nD
  field_names <- x$field_names
  design_names <- x$design_names
  if(x$per_location_design) stop("Not supported: per-location design in `compareGLM`.")
  rm(x)
  # if (verbose>0) {
  #   cat("Number of timepoints:    ",
  #     if (length(unique(nT))==1) { nT } else { paste0(min(nT), "-", max(nT)) }, "\n")
  #   cat("Number of fields:        ", nK, "\n")
  # }

  ### Get `session_names`. -----------------------------------------------------
  session_names <- "single_sess"
  names(BOLD) <- session_names
  names(design) <- session_names

  ### Check `nuisance`. --------------------------------------------------------
  if (!is.null(nuisance)) {
    nuisance <- BayesGLM_format_nuisance(nuisance, nS_expect=nS, nT_expect=nT)

    if (!is.null(names(nuisance)) && !all(names(nuisance) == session_names)) {
      #warning("Ignoring `names(nuisance)`; use `session_names` in `make_design`.")
    }
  } else {
    nuisance <- vector("list", nS)
  }
  names(nuisance) <- session_names

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  valid_cols <- do.call(rbind, lapply(design, function(q) {
    apply(q, 2, function(r){!all(is.na(r))})
  }))
  if (any(colSums(valid_cols)==0)) { stop("Some tasks are missing from every session.") }
  for (ss in seq(nS)) {
    any_bad_design_cols <- if (nD == 1) {
      any(is.na(c(design[[ss]][,valid_cols[ss,]])))
    } else {
      any(is.na(c(design[[ss]][,valid_cols[ss,],])))
    }
    if (any_bad_design_cols) {
      stop("`design` has some sessions & tasks for which some data values ",
        "are `NA`. Partially missing data is not allowed. (Missing tasks ",
        "should have all `NA`.)")
    }
  }

  ### Get `nV`. ----------------------------------------------------------------
  nV <- stop("[TO DO]")

  # QC mask. -------------------------------------------------------------------
  # Mask based on quality control metrics of the BOLD data.
  mask_qc <- make_mask(
    BOLD,
    meanTol=meanTol, varTol=varTol, verbose=verbose>0
  ) #, snrTol=snrTol)
  if (!any(mask_qc$mask)) { stop("No locations meeeting `meanTol` and `varTol`.") }
  if (any(!mask_qc$mask)) {
    if (spatial_type == "mesh") {
      spatial$mask[spatial$mask] <- mask_qc
    } else if (spatial_type == "voxel") {
      spatial$labels[spatial$labels!=0][!mask_qc] <- 0 # [TO DO] check
    } else { stop() }
    BOLD <- lapply(BOLD, function(q){ q[,mask_qc,drop=FALSE] })
  }

  # Update and display the number of data locations. ---------------------------
  nV <- get_nV(spatial, spatial_type)
  if (verbose>0) {
    cat('\tNumber of data locations:', nV$D, '\n')
    if (spatial_type == "voxel") {
      cat('\tNumber of data + boundary locations:', nV$DB, '\n')
    }
  }

  # Scale, nuisance regress, and/or concatenate session data. ------------------

  #collect data and design matrices
  # nK2 <- if (is.null(data[[1]]$nuisance)) { 0 } else { ncol(data[[1]]$nuisance) } #number of nuisance regressors
  for (ss in seq(nS)) {
    # Remove any missing fields from design matrix for classical GLM
    vcols_ss <- valid_cols[ss,]
    # if (!all(vcols_ss)) { warning(
    #   'For session ',ss,', ignoring ',sum(!vcols_ss),
    #   ' design matrix columns of zeros for classical GLM.'
    # )}

    # # Regress nuisance parameters from BOLD data and design matrix.
    # if (!is.null(nuisance[[ss]])) {
    #   nuis_ss <- nuisance[[ss]]
    #   stopifnot((is.matrix(nuis_ss) || is.data.frame(nuis_ss)) && is.numeric(nuis_ss))
    #   nuis_ss <- scale(nuis_ss, scale=FALSE)
    #   # Add intercept to nuisance in case BOLD is not centered.
    #   BOLD[[ss]] <- fMRItools::nuisance_regression(BOLD[[ss]], cbind(1, nuis_ss))
    #   # Do not add intercept, because `design` should already be centered.
    #   design[[ss]][,vcols_ss] <- fMRItools::nuisance_regression(
    #     design[[ss]][,vcols_ss], nuis_ss
    #   ) #[TO DO] if design matrix varies spatially, need to adapt this.
    #   # Design matrix will start as TxKxV and continue in that format after this step.
    #   # [TO DO] Re-scale design?
    #   nuisance[ss] <- list(NULL)
    #   rm(nuis_ss)
    # }

    # Scale data.
    # (`scale_timeseries` expects VxT data, so transpose before and after.)
    BOLD[[ss]] <- t(fMRItools::scale_timeseries(
      t(BOLD[[ss]]), scale=scale_BOLD, transpose=FALSE
    ))
  }
  rm(nuisance)

  # [TO DO] Question: mesh vs surf? same?
  if (spatial_type=="voxel" && do$pw && ar_smooth > 0) {
    stop("[TO DO] How to support this?")
  }
  # mesh$loc,mesh$graph$tv,mesh$n

  stop("TO DO")
  # Replace these:
  # GLM_est_resid_var_pw
  # sparse_and_PW
  # GLM_classical

  x <- GLM_multi(
    BOLD[[1]], design[[1]], nuisance[[1]],
    spatial
  )

  result_classical <- prewhiten_info <- design <- NULL
}
