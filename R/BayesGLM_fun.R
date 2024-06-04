#' BayesGLM_fun
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
#'    \item{nbhd_order}{See documentation for \code{\link{BayesGLM}}.}
#'    \item{buffer}{See documentation for \code{\link{BayesGLM}}.}
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
#' @importFrom fMRItools is_1 nuisance_regression
#'
#' @importFrom utils tail
#'
#' @importFrom methods as
#' @export
BayesGLM_fun <- function(
  BOLD,
  design,
  nuisance=NULL,
  spatial,
  # Below arguments shared with `BayesGLM`.
  scale_BOLD = c("mean", "sd", "none"),
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
  scale_design <- FALSE

  # Initialize return values that may or may not be computed. ------------------
  INLA_model_obj <- hyperpar_posteriors <- Q_theta <- NULL
  field_estimates <- RSS <- hyperpar_posteriors <- mu_theta <- y_all <- XA_all_list <- NULL
  theta_estimates <- Sig_inv <- mesh <- mesh_orig <- NULL

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
  do$Bayesian <- x$Bayes; rm(Bayes) # rename
  do$EM <- x$do_EM; rm(EM) # rename
  do$pw <- x$do_pw # unused
  return_INLA <- x$return_INLA
  rm(x)

  # Modeled after `BayesGLM` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Check `BOLD`. ------------------------------------------------------------
  nS <- length(BOLD)
  ### Check `design`. ----------------------------------------------------------
  # Make `design` a sessions-length list of design matrices.
  #   Get `nK`, `field_names`, and `do$perLocDesign`. Check for consistent dims
  #   across sessions.
  x <- BayesGLM_format_design(design, scale_design=scale_design, nS_expect=nS)
  design <- x$design
  nT <- x$nT #number of fMRI time points
  nK <- x$nK #number of columns in design matrix
  nD <- x$nD
  field_names <- x$field_names
  design_names <- x$design_names
  do$perLocDesign <- x$per_location_design
  rm(x)
  # if (verbose>0) {
  #   cat("Number of timepoints:    ",
  #     if (length(unique(nT))==1) { nT } else { paste0(min(nT), "-", max(nT)) }, "\n")
  #   cat("Number of fields:        ", nK, "\n")
  # }

  ### Get `session_names`. -----------------------------------------------------
  session_names <- BayesGLM_session_names(
    nS, session_names=NULL, names(BOLD), names(design)
  )
  names(BOLD) <- session_names
  names(design) <- session_names

  # if (verbose>0) {
  #   cat("Session names:           ", paste0(session_names, collapse=", "), "\n")
  #   cat("Field names:             ", paste0(field_names, collapse=", "), "\n")
  # }

  ### Check `nuisance`. --------------------------------------------------------
  if (!is.null(nuisance) && !all(vapply(nuisance, is.null, FALSE))) {
    nuisance <- BayesGLM_format_nuisance(nuisance, nS_expect=nS, nT_expect=nT)

    if (!is.null(names(nuisance)) && !all(names(nuisance) == session_names)) {
      #warning("Ignoring `names(nuisance)`; use `session_names` in `make_design`.")
    }
  } else {
    nuisance <- vector("list", nS)
  }
  names(nuisance) <- session_names

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  design_type <- if (do$perLocDesign) { "per_location" } else { "regular" }
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

  ### Check `spatial` and get `nV`. --------------------------------------------

  stopifnot(is.list(spatial))
  mesh_spatial_names <- c("surf", "mask")
  is_mesh <- length(names(spatial))==2 && all(names(spatial)==mesh_spatial_names)
  voxel_spatial_names <- c(
    "labels", "trans_mat", "trans_units", "nbhd_order", "buffer", "buffer_mask"
  )
  is_voxel <- length(names(spatial))==6 && all(names(spatial)==voxel_spatial_names)
  if (!is_mesh && !is_voxel) { stop("`spatial` is not correctly formatted. Please fix.") }
  spatial_type <- if (is_mesh) { "mesh" } else { "voxel" }
  rm(is_mesh, is_voxel, mesh_spatial_names, voxel_spatial_names)
  nV <- get_nV(spatial, spatial_type) #total number of locations, number of data locations

  if (spatial_type=="mesh" && (nV$T != nV$D)) {
    if (verbose>0) {
      cat(paste0(
        "\t", (nV$T-nV$D),
        " locations on the mesh do not have data.\n"
      ))
    }
  }

  # QC mask. -------------------------------------------------------------------
  # Mask based on quality control metrics of the BOLD data.
  mask_qc <- make_mask(
    BOLD,
    meanTol=meanTol, varTol=varTol, verbose=verbose>0
  ) #, snrTol=snrTol)
  if (!any(mask_qc$mask)) { stop("No locations meeeting `meanTol` and `varTol`.") }
  if (any(!mask_qc$mask)) {
    if (spatial_type == "mesh") {
      spatial$mask[spatial$mask] <- mask_qc$mask
    } else if (spatial_type == "voxel") {
      spatial$labels[spatial$labels!=0][!mask_qc$mask] <- 0
    } else { stop() }
    BOLD <- lapply(BOLD, function(q){ q[,mask_qc$mask,drop=FALSE] })
  }

  # Get SPDE and mask based on it (additional vertices may be excluded).
  x <- switch(spatial_type, mesh=SPDE_from_mesh, voxel=SPDE_from_voxel)(spatial)
  if (spatial_type=="mesh") {
    BOLD <- lapply(BOLD, function(q){ q[,x$mask_new_diff,drop=FALSE] })
  }
  spde <- x$spde
  spatial <- x$spatial
  rm(x)

  # Adjust design for per-location modeling.

  # Update and display the number of data locations. ---------------------------
  # (Need to update after QC mask and SPDE)
  nV <- get_nV(spatial, spatial_type)
  if (verbose>0) {
    cat('\tNumber of data locations:', nV$D, '\n')
    if (spatial_type == "voxel") {
      cat('\tNumber of data + boundary locations:', nV$DB, '\n')
    }
  }

  # Scale, nuisance regress, and/or concatenate session data. ------------------

  # Check for intercepts and flat columns in design and nuisance matrices.
  des_has_intercept <- nuis_has_intercept <- vector("logical", nS)
  for (ss in seq(nS)) {
    # Stop if any zero-var, zero-mean column exists.
    des_ss_is_flat <- apply(abs(design[[ss]]) < 1e-8, 2, all)
    if (any(des_ss_is_flat)) {
      stop("Design matrix for session ", ss, " has at least one column that is ",
        "flat (all values are near-zero).")
    }

    # Detect zero-var, nonzero-mean columns.
    if (design_type=="per_location") {
      des_ss_is_intercept <- apply(design[[ss]], c(2,3), var) < 1e-8
    } else {
      des_ss_is_intercept <- matrixStats::colVars(design[[ss]]) < 1e-8
    }
    if (sum(des_ss_is_intercept) > 1) {
      stop("Design matrix for session ", ss, " has more than one intercept-like ",
        "column, which will cause problems during model fitting. Please fix.")
    }
    des_has_intercept[ss] <- any(des_ss_is_intercept)

    # For nuisance: detect zero-var, nonzero-mean columns.
    if (!is.null(nuisance[[ss]])) {
      nuis_ss_is_intercept <- matrixStats::colVars(nuisance[[ss]]) < 1e-8
      nuis_has_intercept[ss] <- any(nuis_ss_is_intercept)
    } else {
      nuis_has_intercept[ss] <- FALSE
    }
  }
  if(any(des_has_intercept)) stop('Design matrix must not have an intercept, since data and design will be centered prior to model fitting. Please fix.')

  #collect data and design matrices
  nK2 <- vector("numeric", nS)
  for (ss in seq(nS)) {
    # Remove any missing fields from design matrix for classical GLM
    vcols_ss <- valid_cols[ss,]
    # if (!all(vcols_ss)) { warning(
    #   'For session ',ss,', ignoring ',sum(!vcols_ss),
    #   ' design matrix columns of zeros for classical GLM.'
    # )}

    # THIS PART HAS BEEN MODIFIED BECAUSE DESIGN NO LONGER ALLOWED TO HAVE INTERCEPT
    # If the design matrix has an intercept, ensure the nuisance matrix does not.
    # If neither have an intercept, add one to the nuisance for centering
    #   `BOLD` and `design`.
    #if (des_has_intercept[ss]) {
      # if (nuis_has_intercept[ss]) {
      #   stop("Both the design and nuisance for session ", ss, " has an ",
      #     "intercept column. Remove one.")
      # }
      # if (verbose > 0) {
      #   cat(paste0("\t\tSession ", ss, ": modeling BOLD means using the intercept field in the provided `design`.\n"))
      # }
    #} else
    if (!nuis_has_intercept[ss]) {
      nuisance[[ss]] <- cbind2(nuisance[[ss]], as.matrix(rep(1, nT[ss])))
      nuis_has_intercept[ss] <- TRUE
    }
    if (verbose > 0) { if(ss==1) cat("\tDemeaning BOLD and design during nuisance regression.\n") }

    # NOTE: At this point, it is not possible for nuisance[[ss]] to be NULL.

    # Regress nuisance parameters from BOLD data and design matrix.
    nK2[ss] <- if (is.null(nuisance[[ss]])) { 0 } else { ncol(nuisance[[ss]]) }
    if (!is.null(nuisance[[ss]])) {
      nuis_ss <- nuisance[[ss]]
      stopifnot((is.matrix(nuis_ss) || is.data.frame(nuis_ss)) && is.numeric(nuis_ss))
      if (!nuis_has_intercept[ss]) {
        stop('nuisance matrix has no intercept, something is wrong. Contact developer.')
        #nuis_ss <- scale(nuis_ss, scale=FALSE)
      }

      #center, filter/detrend, and nuisance regress the BOLD
      BOLD_mean_ss <- apply(BOLD[[ss]], MARGIN = 2, FUN = median) # [TO DO] instead, save the intercept from nuisance regression
      BOLD[[ss]] <- fMRItools::nuisance_regression(BOLD[[ss]], nuis_ss)
      #center, filter/detrend, and nuisance regress the BOLD
      if (!do$perLocDesign) {
        design[[ss]][,vcols_ss] <- fMRItools::nuisance_regression(
          design[[ss]][,vcols_ss], nuis_ss
        )
      } else {
        design[[ss]][,,vcols_ss][] <- fMRItools::nuisance_regression(
          matrix(design[[ss]][,,vcols_ss], nrow=dim(design[[ss]])[1]),
          nuis_ss
        )
      }
      nuisance[ss] <- list(NULL)
      rm(nuis_ss)
    } else {
      stop('nuisance matrix is NULL, something is wrong. Contact developer.')
    }

    # Scale data.
    # (`scale_BOLD` expects VxT data, so transpose before and after.)
    BOLD[[ss]] <- t(scale_BOLD(t(BOLD[[ss]]), scale=scale_BOLD, v_means = BOLD_mean_ss))

  }
  rm(des_has_intercept, nuis_has_intercept, vcols_ss, nuisance)

  # Estimate residual variance (for var. std.izing) and get prewhitening info.
  if (do$pw && verbose>0) { cat("\tEstimating prewhitening parameters.\n") }
  x <- GLM_est_resid_var_pw(
    BOLD, design, spatial, spatial_type,
    session_names, field_names, design_type,
    valid_cols, nT, nD,
    ar_order, ar_smooth, aic, n_threads,
    do$pw
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
      BOLD[[ss]], design[[ss]],
      spatial, spatial_type,
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
      BOLD[[ss]], design[[ss]], nK2[ss], nV$D,
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
      #post-multiply each design matrix by A (n_data x nMesh) for Bayesian GLM
      XA_ss <- lapply(design[[ss]], function(x) { x %*% A_sparse_ss })
      #cbind (expanded) design matrices for each field
      XA_ss <- list(do.call(cbind, XA_ss))
      names(XA_ss) <- session_names[ss]
      XA_all_list <- c(XA_all_list, XA_ss)
      #rm(XA_ss, A_sparse_ss)
    }
  }

  # [NOTE] Moved to `GLM_FIR.R`: FIR Model.

  # Bayesian GLM. --------------------------------------------------------------
  if (do$Bayesian) {

    # Construct betas and repls objects.
    x <- make_replicates(
      nSess=nS, field_names=field_names, nMesh=nV$D
      #, data_loc=data_loc) #indices of original data locations
    )
    betas <- x$betas
    repls <- x$repls
    rm(x)

    model_data <- make_data_list(y=y_all, X=XA_all_list, betas=betas, repls=repls)
    Amat <- model_data$X
    model_data$XA_all_list <- NULL

    # [NOTE] Moved to `GLM_Bayesian_EM.R`: EM Model.

    ## INLA Model. -------------------------------------------------------------
    #estimate model using INLA
    if (verbose>0) cat('\tEstimating Bayesian model with INLA...')
    #organize the formula and data objects
    hyper_initial <- c(-2,2)
    hyper_initial <- rep(list(hyper_initial), nK)
    hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

    formula <- paste0('f(',field_names, ', model = spde, replicate = ', names(repls), hyper_vec, ')')
    formula <- paste(c('y ~ -1', formula), collapse=' + ')
    formula <- as.formula(formula)

    INLA_model_obj <- INLA::inla(
      formula,
      data=model_data,
      #data=INLA::inla.stack.data(model_data, spde=spde),
      control.predictor=list(A=Amat, compute = TRUE),
      verbose = verbose>1, keep = FALSE, num.threads = n_threads,
      control.inla = list(strategy = "gaussian", int.strategy = "eb"),
      control.family=list(hyper=list(prec=list(initial=1))),
      control.compute=list(config=TRUE), contrasts = NULL, lincomb = NULL #required for excursions
    )
    if (verbose>0) cat("\tDone!\n")

    #extract useful stuff from INLA model result
    field_estimates <- extract_estimates(
      INLA_model_obj=INLA_model_obj,
      session_names=session_names,
      spatial=spatial, spatial_type=spatial_type, spde=spde
    ) #posterior means of latent field
    theta_estimates <- INLA_model_obj$summary.hyperpar$mean
    RSS <- vector("list", nS)
    for (ss in seq(nS)) {
      RSS[[ss]] <- rowSums(t(matrix(
        as.matrix(c(BOLD[[ss]])) - (
          do.call(cbind, design[[ss]][valid_cols[ss,]]) %*% as.matrix(c(field_estimates[[ss]]))
        ),
        nrow = nT[ss]
      ))^2)
    }
    # resids <- t(matrix(as.matrix(y) - (X %*% coefs), nrow = nT_ss))
    hyperpar_posteriors <- get_posterior_densities2(
      INLA_model_obj=INLA_model_obj, #spde,
      field_names
    ) #hyperparameter posterior densities

    #construct object to be returned
    INLA_model_obj <- switch(return_INLA,
      trimmed=trim_INLA_model_obj(INLA_model_obj, minimal=FALSE),
      full=INLA_model_obj,
      minimal=trim_INLA_model_obj(INLA_model_obj, minimal=TRUE)
    )
    attr(INLA_model_obj, "format") <- return_INLA
  } else {
    field_estimates <- lapply(result_classical, function(x){ x$estimates })
    attr(field_estimates, "GLM_type") <- "classical"
  }

  result <- list(
    field_estimates = field_estimates, # new in 6.0: these are masked.
    INLA_model_obj = INLA_model_obj,
    hyperpar_posteriors = hyperpar_posteriors,
    theta_estimates = theta_estimates,
    RSS = RSS,
    result_classical = result_classical,
    #result_FIR = result_FIR,
    spatial = spatial,
    spde = spde,
    field_names = field_names,
    session_names = session_names,
    mask_qc = mask_qc,
    # For joint group model ~~~~~~~~~~~~~
    #posterior_Sig_inv = Sig_inv,
    y = y_all,
    X = XA_all_list,
    prewhiten_info = prewhiten_info,
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call = match.call()
  )
  class(result) <- "BGLM0"

  result
}
