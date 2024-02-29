#' BayesGLM
#'
#' Performs spatial Bayesian GLM for task fMRI activation
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
#' @return A \code{"BayesGLM"} object: a list with elements
#'  \describe{
#'    \item{INLA_model_obj}{The full result of the call to \code{INLA::inla}.}
#'    \item{field_estimates}{The estimated coefficients for the Bayesian model.}
#'    \item{result_classical}{Results from the classical model: field estimates, field standard error estimates, residuals, degrees of freedom, and the mask.}
#'    \item{mesh}{The model mesh including only the locations analyzed, i.e. within \code{mask}, without missing values, and meeting \code{meanTol} and \code{varTol}.}
#'    \item{mesh_orig}{The original mesh provided.}
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
BayesGLM <- function(
  BOLD, design, nuisance=NULL,
  spatial,
  # Below arguments shared with `BayesGLM_cifti`.
  scale_BOLD = c("auto", "mean", "sd", "none"),
  Bayes = TRUE,
  #EM = FALSE,
  ar_order = 6,
  ar_smooth = 5,
  aic = FALSE,
  n_threads = 4,
  return_INLA = c("trimmed", "full", "minimal"),
  verbose = 1,
  meanTol = 1e-6,
  varTol = 1e-6#,
  #snrTol = 50,
  #emTol = 1e-3
  ){

  EM <- FALSE
  emTol <- 1e-3

  # Initialize return values that may or may not be computed. ------------------
  INLA_model_obj <- hyperpar_posteriors <- Q_theta <- NULL
  field_estimates <- hyperpar_posteriors <- mu_theta <- y_all <- XA_all_list <- NULL
  theta_estimates <- Sig_inv <- mesh <- mesh_orig <- NULL

  # Preliminaries. -------------------------------------------------------------
  ### Simple parameter checks. -------------------------------------------------

  session_names <- names(design$design)
  nS <- nS_orig <- length(session_names)
  field_names <- design$field_names
  nK <- length(field_names)

  design_type <- design$design_type
  valid_cols <- design$valid_cols
  design <- design$design

  nT <- vapply(design, nrow, 0)
  nD <- if (design_type %in% c("design", "onsets")) { 1 } else { dim(design[[1]])[3] }

  stopifnot(is.list(spatial))
  mesh_spatial_names <- c("surf", "mask")
  is_mesh <- length(names(spatial))==2 && all(names(spatial)==mesh_spatial_names)
  voxel_spatial_names <- c(
    "label", "trans_mat", "trans_units", "nbhd_order", "buffer"
  )
  is_voxel <- length(names(spatial))==5 && all(names(spatial)==voxel_spatial_names)
  if (!is_mesh && !is_voxel) { stop("`spatial` is not correctly formatted. Please fix.") }
  spatial_type <- if (is_mesh) { "mesh" } else { "voxel" }
  rm(is_mesh, is_voxel, mesh_spatial_names, voxel_spatial_names)
  nV <- get_nV(spatial, spatial_type)

  do <- vector("list")

  # In a separate function because these checks are shared with `BayesGLM_cifti`.
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
  do$Bayesian <- x$Bayes; rm(Bayes) # rename
  do$EM <- x$do_EM; rm(EM) # rename
  do$pw <- x$do_pw
  return_INLA <- x$return_INLA
  rm(x)

  # More checks, in case `BayesGLM` was run directly (not w/ CIFTI) ------------
  # [TO DO] Refine.
  if (is.null(nuisance)) { nuisance <- vector("list", nS) }
  for (ss in seq(nS)) {
    stopifnot(is.matrix(BOLD[[ss]]) && is.numeric(BOLD[[ss]]))
    stopifnot(nrow(BOLD[[ss]]) == nT[ss])
    stopifnot(ncol(BOLD[[ss]]) == nV$D)
    if (!is.null(nuisance[[ss]])) {
      stopifnot(nrow(nuisance[[ss]]) == nT[ss])
    }
  }

  ## Set up spatial data, SPDE -------------------------------------------------
  spatial_og <- spatial

  # Mask based on quality control metrics of the BOLD data.
  mask_qc <- make_mask(BOLD, meanTol=meanTol, varTol=varTol) #, snrTol=snrTol)
  if (!any(mask_qc$mask)) { stop("No locations meeeting `meanTol` and `varTol`.") }
  if (any(!mask_qc$mask)) {
    if (spatial_type == "mesh") {
      spatial$mask[spatial$mask] <- mask_qc
    } else if (spatial_type == "voxel") {
      spatial$labels[spatial$labels!=0][!mask_qc] <- 0 # [TO DO] check
    } else { stop() }
    BOLD <- lapply(BOLD, function(q){ q[,mask_qc,drop=FALSE] })
  }

  # Get SPDE and mask based on it (additional vertices may be excluded).
  x <- switch(spatial_type, mesh=SPDE_from_mesh, voxel=SPDE_from_voxel)(spatial)
  if (spatial_type=="mesh") {
    BOLD <- lapply(BOLD, function(q){ q[,x$mask_new_diff,drop=FALSE] })
  }
  spde <- x$spde
  spatial <- x$spatial
  rm(x)

  # Get and display the number of data locations. ------------------------------
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
    if (!all(vcols_ss)) { warning(
      'For session ',ss,', ignoring ',sum(!vcols_ss),
      ' design matrix columns of zeros for classical GLM.'
    )}

    # Regress nuisance parameters from BOLD data and design matrix.
    if (!is.null(nuisance[[ss]])) {
      nuis_ss <- nuisance[[ss]]
      stopifnot((is.matrix(nuis_ss) || is.data.frame(nuis_ss)) && is.numeric(nuis_ss))
      nuis_ss <- scale(nuis_ss, scale=FALSE)
      # Add intercept to nuisance in case BOLD is not centered.
      BOLD[[ss]] <- fMRItools::nuisance_regression(BOLD[[ss]], cbind(1, nuis_ss))
      # Do not add intercept, because `design` should already be centered.
      design[[ss]][,vcols_ss] <- fMRItools::nuisance_regression(
        design[[ss]][,vcols_ss], nuis_ss
      ) #[TO DO] if design matrix varies spatially, need to adapt this.
      # Design matrix will start as TxKxV and continue in that format after this step.
      # [TO DO] Re-scale design?
      nuisance[ss] <- list(NULL)
    }

    # Scale data.
    # (`scale_timeseries` expects VxT data, so transpose before and after.)
    BOLD[[ss]] <- t(fMRItools::scale_timeseries(
      t(BOLD[[ss]]), scale=scale_BOLD, transpose=FALSE
    ))
  }
  rm(nuisance, nuis_ss)

  # [TO DO] Question: mesh vs surf? same?
  if (spatial_type=="voxel" && do$pw && ar_smooth > 0) {
    stop("[TO DO] How to support this?")
  }
  # mesh$loc,mesh$graph$tv,mesh$n

  # Estimate residual variance (for var. std.izing) and get prewhitening info.
  x <- GLM_est_resid_var_pw(
    BOLD, design, spatial, spatial_type,
    session_names, field_names, design_type,
    valid_cols, nT, nD,
    ar_order, ar_smooth, aic, n_threads,
    do$pw, verbose
  )
  var_resid <- x$var_resid
  sqrtInv_all <- x$sqrtInv_all # `NULL` if `!do$pw`
  prewhiten_info <- x[c("AR_coefs_avg", "var_avg", "max_AIC")]
  rm(x)

  # Classical GLM. -------------------------------------------------------------

  result_classical <- setNames(vector('list', length=nS), session_names)
  for (ss in seq(nS)) {
    if (verbose>0) {
      if (nS==1) {
        cat('\tFitting classical GLM.\n')
      } else {
        if (ss==1) { cat('\tFitting classical GLM for session ') }
        if (ss!=nS) { cat(paste0(ss, ", ")) } else { cat(paste0(ss, ".\n")) }
      }
    }

    vcols_ss <- valid_cols[ss,]

    # Set up vectorized data and big sparse design matrix.
    # Apply prewhitening, if applicable.
    x <- sparse_and_PW(
      BOLD[[ss]], design[[ss]], nV$T, nV$D,
      field_names, design_type,
      vcols_ss, nT[ss], nD,
      sqrtInv_all[[ss]]
    )
    BOLD[[ss]] <- x$BOLD
    design[[ss]] <- x$design
    A_sparse_ss <- x$A_sparse
    rm(x)

    # Compute classical GLM.
    result_classical[[ss]] <- GLM_classical(
      BOLD[[ss]], design[[ss]], nV$D,
      field_names, design_type,
      vcols_ss, nT[ss], nD,
      var_resid, sqrtInv_all[[ss]],
      do$pw, compute_SE=TRUE
    )

    # Set up for Bayesian GLM.
    if (ss==1) {
      y_all <- vector("numeric")
      XA_all_list <- NULL
    }
    if (do$Bayesian) {
      y_all <- c(y_all, BOLD[[ss]])
      #XA_ss <- design
      #post-multiply each design matrix by A (n_data x n_mesh) for Bayesian GLM
      XA_ss <- lapply(design, function(x) { x %*% A_sparse_ss })
      #cbind (expanded) design matrices for each field
      XA_ss <- list(do.call(cbind, XA_ss))
      names(XA_ss) <- session_names[ss]
      XA_all_list <- c(XA_all_list, XA_ss)
      rm(XA_ss, A_sparse_ss)
    }
  }

  # [NOTE] Moved to `GLM_FIR.R`: FIR Model.

  # # Bayesian GLM. --------------------------------------------------------------
  # if (do$Bayesian) {
  #   #construct betas and repls objects
  #   if (spatial_type == "mesh") {
  #     n_mesh <- length(mesh$idx$loc)
  #     if(!all.equal(mesh$idx$loc, seq(n_mesh))) stop('Developer: Check definition of `spatial` in organize_replicates')
  #     #data_loc <- seq(n_mesh)
  #   } else if (spatial_type == "voxel") {
  #     n_mesh <- spde$n.spde
  #     #data_loc <- data_loc
  #   } else { stop() }

  #   x <- organize_replicates(
  #     n_sess=nS, field_names=field_names, n_mesh=n_mesh#, data_loc=data_loc) #indices of original data locations
  #   )
  #   betas <- x$betas
  #   repls <- x$repls
  #   rm(x)

  #   model_data <- make_data_list(y=y_all, X=XA_all_list, betas=betas, repls=repls)
  #   Amat <- model_data$X
  #   model_data$XA_all_list <- NULL

  #   # [NOTE] Moved to `GLM_Bayesian_EM.R`: Em Model.

  #   ## INLA Model. -------------------------------------------------------------
  #   #estimate model using INLA
  #   if (verbose>0) cat('\tEstimating Bayesian model with INLA...')
  #   #organize the formula and data objects
  #   hyper_initial <- c(-2,2)
  #   hyper_initial <- rep(list(hyper_initial), nK)
  #   hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

  #   formula <- paste0('f(',field_names, ', model = spde, replicate = ', names(repls), hyper_vec, ')')
  #   formula <- paste(c('y ~ -1', formula), collapse=' + ')
  #   formula <- as.formula(formula)

  #   INLA_model_obj <- INLA::inla(
  #     formula,
  #     data=model_data,
  #     #data=INLA::inla.stack.data(model_data, spde=spde),
  #     control.predictor=list(A=Amat, compute = TRUE),
  #     verbose = verbose>1, keep = FALSE, num.threads = n_threads,
  #     control.inla = list(strategy = "gaussian", int.strategy = "eb"),
  #     control.family=list(hyper=list(prec=list(initial=1))),
  #     control.compute=list(config=TRUE), contrasts = NULL, lincomb = NULL #required for excursions
  #   )
  #   if (verbose>0) cat("\tDone!\n")

  #   #extract useful stuff from INLA model result
  #   field_estimates <- extract_estimates(
  #     INLA_model_obj=INLA_model_obj,
  #     session_names=session_names,
  #     mask=nV$maskD, # [TO DO]
  #     inds=nV # [TO DO]
  #   ) #posterior means of latent field
  #   theta_estimates <- INLA_model_obj$summary.hyperpar$mean
  #   hyperpar_posteriors <- get_posterior_densities2(
  #     INLA_model_obj=INLA_model_obj, #spde,
  #     field_names
  #   ) #hyperparameter posterior densities

  #   #construct object to be returned
  #   INLA_model_obj <- switch(return_INLA,
  #     trimmed=trim_INLA_model_obj(INLA_model_obj, minimal=FALSE),
  #     full=INLA_model_obj,
  #     minimal=trim_INLA_model_obj(INLA_model_obj, minimal=TRUE)
  #   )
  #   attr(INLA_model_obj, "format") <- return_INLA
  # } else {
    field_estimates <- lapply(result_classical, function(x){ x$estimates })
    attr(field_estimates, "GLM_type") <- "classical"
  # }

  result <- list(
    field_estimates = field_estimates, # new in 6.0: these are masked.
    INLA_model_obj = INLA_model_obj,
    hyperpar_posteriors = hyperpar_posteriors,
    theta_estimates = theta_estimates,
    result_classical = result_classical,
    #result_FIR = result_FIR,
    spatial = spatial,
    spde = spde,
    mask_qc = mask_qc,
    field_names = field_names,
    session_names = session_names,
    # For joint group model ~~~~~~~~~~~~~
    #posterior_Sig_inv = Sig_inv,
    y = y_all,
    X = XA_all_list,
    prewhiten_info = prewhiten_info,
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call = match.call()
  )
  class(result) <- "BayesGLM"

  result
}
