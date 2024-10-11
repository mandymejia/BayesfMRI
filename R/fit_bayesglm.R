#' fit_bayesglm
#'
#' Performs spatial Bayesian GLM for task fMRI activation
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param BOLD,design,nuisance Session-length list of numeric matrices/arrays,
#'  each with volumes along the first dimension.
#' @param scrub Session-length list of spike regressors: numeric matrices, with
#'  volumes along the first dimension, valued at 1 for scrubbed volumes and 0
#'  otherwise.
#'
#'  Scrubbing is performed by incorporating spike regressors in the nuisance
#'  matrix during nuisance regression (in a simultaneous framework), and then
#'  removing the scrubbed timepoints from the resulting BOLD and design.
#' @param spatial Gives the spatial information:
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
#' @param hyperpriors Should informative or default non-informative hyperpriors be assumed on SPDE hyperparameters?
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
#'    \item{mesh}{The model mesh.}
#'    \item{mask}{A mask of \code{mesh} indicating the locations inside \code{mesh}.}
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
#' @importFrom stats as.formula var dist
#' @importFrom fMRItools is_1 nuisance_regression
#'
#' @importFrom utils tail
#'
#' @importFrom methods as
#' @export
fit_bayesglm <- function(
  BOLD,
  design,
  nuisance=NULL,
  scrub=NULL,
  spatial,
  # Below arguments shared with `BayesGLM`.
  scale_BOLD = c("mean", "sd", "none"),
  Bayes = TRUE,
  hyperpriors = c("informative","default"),
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
  hyperpriors <- hyperpriors[1]
  if(!hyperpriors %in% c("informative","default")) stop('`hyperpriors` must be "informative" or "default"')

  # Initialize return values that may or may not be computed. ------------------
  INLA_model_obj <- hyperpar_posteriors <- Q_theta <- NULL
  field_ests <- RSS <- hyperpar_posteriors <- mu_theta <- y_all <- XA_all_list <- NULL
  theta_ests <- theta_ests2 <- Sig_inv <- mesh <- NULL

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
  # nT <- vapply(design, nrow, 0) #number of fMRI time points # not used; wait till after scrubbing
  nK <- dim(design[[1]])[2] #number of columns in design matrix
  field_names <- colnames(design[[1]]) #names of task regressors
  do$perLocDesign <- (length(dim(design[[1]])) == 3)
  #nD <- if(do$perLocDesign) sapply(design, function(x) dim(x)[3]) else rep(1, nS) #number of design matrices, if allowed to vary by location
  session_names <- names(design)

  ### Validate `spatial`. ------------------------------------------------------
  validate_spatial(spatial)
  spatial_type <- spatial$spatial_type

  # Initialize and construct informative priors for SPDE hyperparameters -------
  logkappa_vec <- logtau_vec <- logkappa0 <- logtau0 <- NULL # default priors
  if (do$Bayesian) {
    x <- log_kappa_tau(spatial, hyperpriors, verbose)
    logkappa_vec <- x$logkappa_vec
    logkappa0 <- x$logkappa0
    logtau <- x$logtau
    logtau0 <- x$logtau0
  }

  # QC mask, priors, mesh, and SPDE. -------------------------------------------

  # Print message about medial wall mask, for cortex model.
  nV <- get_nV(spatial)
  if (spatial_type=="vertex" && (nV$total != nV$input)) {
    if (verbose>0) {
      cat(paste0("\t", (nV$total-nV$input), " vertices do not have data.\n"))
    }
  }

  # Get QC mask.
  BOLD_QC <- do_QC(
    BOLD, meanTol=meanTol, varTol=varTol, verbose=verbose>0
  ) #, snrTol=snrTol)
  if (!any(BOLD_QC$mask)) { stop("No locations meeeting `meanTol` and `varTol`.") }

  # Apply QC mask to BOLD data.
  if (any(!BOLD_QC$mask)) {
    BOLD <- lapply(BOLD, function(q){ q[,BOLD_QC$mask,drop=FALSE] }) #remove bad locations from BOLD
  }

  # Update `spatial` w/ QC info and mesh for cortex model.
  # Get SPDE, using `spatial` and the QC mask.
  x <- switch(spatial_type, vertex=SPDE_from_vertex, voxel=SPDE_from_voxel)(
    spatial, qc_mask = BOLD_QC$mask,
    logkappa = logkappa_vec, logtau = logtau_vec
  )
  spde <- x$spde
  spatial <- x$spatial
  rm(x)

  # Update and display the number of data locations. ---------------------------
  nV <- get_nV(spatial)
  if (verbose>0) {
    cat('\tNumber of modeled data locations:', nV$mdata, '\n')
    cat('\tNumber of modeled data + boundary locations:', nV$model, '\n')
  }

  # Nuisance regression, scrubbing, and scaling. -------------------------------
  ### Do nuisance regression of BOLD and design, and BOLD scaling. -------------

  #identify valid design matrix columns for each session
  valid_cols <- do.call(rbind, lapply(design, function(q) {
    apply(q, 2, function(r){!all(is.na(r))})
  }))

  nK2 <- vector("numeric", nS)
  scrub_vec <- vector("list", nS)
  for (ss in seq(nS)) {
    if (verbose > 0) { if(ss==1) cat("\tBOLD and design nuisance regression.\n") }

    # Remove any missing fields from design matrix for classical GLM
    vcols_ss <- valid_cols[ss,]

    # Merge scrubbing spike regressors with the nuisance matrix
    nuis_ss <- cbind2(scrub[[ss]], nuisance[[ss]])
    # nK2[ss] <- ncol(nuis_ss)
    stopifnot((is.matrix(nuis_ss) || is.data.frame(nuis_ss)) && is.numeric(nuis_ss))

    # Regress nuisance parameters from BOLD data and design matrix.
    #center, filter/detrend, and nuisance regress the BOLD
    BOLD_mean_ss <- apply(BOLD[[ss]], MARGIN = 2, FUN = median) # [TO DO] instead, save the intercept from nuisance regression
    BOLD[[ss]] <- fMRItools::nuisance_regression(BOLD[[ss]], nuis_ss)
    #center, filter/detrend, and nuisance regress the design
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

    # Scrub BOLD and design.
    if (!is.null(scrub[[ss]])) {
      # Get the timepoints to remove.
      scrub_vec[[ss]] <- apply(scrub[[ss]] != 0, 1, any)

      # Remove. # [TO DO] fix
      cat(paste0("\tScrubbing ", sum(scrub_vec[[ss]]),
        " volumes from session ", ss, ".\n"))
      BOLD[[ss]] <- BOLD[[ss]][!scrub_vec[[ss]],,drop=FALSE]
      design[[ss]] <- design[[ss]][!scrub_vec[[ss]],,drop=FALSE]
    }

    # Scale data.
    # (`scale_BOLD` expects VxT data, so transpose before and after.)
    # [TO DO] Scale based off scrubbed data, but scale all data points.
    BOLD[[ss]] <- t(
      scale_BOLD(t(BOLD[[ss]]),
      scale=scale_BOLD, v_means = BOLD_mean_ss)
    )
  }
  rm(vcols_ss, nuisance)

  # Update `nT` due to scrubbing.
  nT <- vapply(BOLD, nrow, 0)

  # Estimate residual variance (for var. std.izing) and get prewhitening info.
  if (do$pw && verbose > 0) {
    x <- if (ar_smooth > 0) { " and smoothing" } else { "" }
    cat(paste0("\tEstimating", x, " prewhitening parameters."))
  }

  design_type <- if (do$perLocDesign) { "per_location" } else { "regular" }
  x <- GLM_est_resid_var_pw(
    BOLD, design, spatial,
    session_names, field_names, design_type,
    valid_cols, nT,
    ar_order, ar_smooth, aic, n_threads, do$pw
  )
  var_resid <- x$var_resid
  sqrtInv_all <- x$sqrtInv_all # diagonal if !do$pw
  prewhiten_info <- x[c("AR_coefs_avg", "var_avg", "max_AIC", "sqrtInv_all")]
  rm(x)

  # Classical GLM. -------------------------------------------------------------
  result_classical <- setNames(vector('list', length=nS), session_names)
  for (ss in seq(nS)) {
    if (verbose>0) {
      if (nS==1) {
        cat('\n\tFitting classical GLM.\n')
      } else {
        if (ss==1) { cat('\n\tFitting classical GLM for session ') }
        if (ss!=nS) { cat(paste0(ss, ", ")) } else { cat(paste0(ss, ".\n")) }
      }
    }

    vcols_ss <- valid_cols[ss,]

    # Set up vectorized data and big sparse design matrix.
    # Apply prewhitening, if applicable.
    x <- sparse_and_PW(
      BOLD[[ss]], design[[ss]],
      spatial, spde,
      field_names, design_type,
      vcols_ss, nT[ss],
      sqrtInv_all[[ss]]
    )
    BOLD[[ss]] <- x$BOLD
    design[[ss]] <- x$design
    A_sparse_ss <- x$A_sparse
    rm(x)

    # Compute classical GLM.
    result_classical[[ss]] <- GLM_classical(
      BOLD[[ss]], design[[ss]], nK2[ss], nV$mdata,
      field_names, design_type,
      vcols_ss, nT[ss],
      do$pw, compute_SE=TRUE
    )

    # #disabled this because it is very close to 1 after prewhitening
    # s2_init <- mean(apply(result_classical[[ss]]$resids, 1, var), na.rm=TRUE)
    # print(paste0('initial value for precision: 1 / ', s2_init))

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

  # Bayesian GLM. --------------------------------------------------------------
  if (do$Bayesian) {

    # Construct betas and repls objects.
    x <- make_replicates(nSess=nS, field_names=field_names, spatial)
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
    #hyper_initial <- c(-2,2)
    #hyper_initial <- rep(list(hyper_initial), nK)
    hyper_initial <- round(c(logkappa0, logtau0),2)
    hyper_initial <- rep(list(hyper_initial), nK)
    hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

    formula <- paste0('f(',field_names, ', model = spde, replicate = ', names(repls), hyper_vec, ')')
    #formula <- paste0('f(',field_names, ', model = spde, replicate = ', names(repls), ')')
    formula <- paste(c('y ~ -1', formula), collapse=' + ')
    formula <- as.formula(formula)

    INLA_model_obj <- INLA::inla(
      formula,
      data=model_data,
      #data=INLA::inla.stack.data(model_data, spde=spde),
      control.predictor=list(A=Amat, compute = TRUE),
      verbose = verbose>1, keep = FALSE, num.threads = n_threads,
      control.inla = list(strategy = "gaussian", int.strategy = "eb"),
      control.family=list(
        hyper=list(prec=list(initial=0, #log(1) = 0
        param=c(1, 1)))
      ), #put a more informative prior on the residual precision
      control.compute=list(config=TRUE), contrasts = NULL, lincomb = NULL #required for excursions
    )
    if (verbose>0) cat("\tDone!\n")

    # Extract stuff from INLA model result -------------------------------------
    field_ests <- extract_estimates(
      INLA_model_obj=INLA_model_obj,
      session_names=session_names,
      spatial=spatial, spde=spde
    ) #posterior means of latent field

    theta_ests <- INLA_model_obj$summary.hyperpar$mode
    names(theta_ests) <- rownames(INLA_model_obj$summary.hyperpar)
    names(theta_ests) <- gsub("Theta1", "log_tau",  names(theta_ests))
    names(theta_ests) <- gsub("Theta2", "log_kappa",  names(theta_ests))

    #compute spatial range and var for each task
    d <- if (spatial_type == "vertex"){ 2 } else { 3 }
    nu <- 2 - d/2
    theta_ests2 <- matrix(NA, length(field_names), 2)
    colnames(theta_ests2) <- c('range','var')
    rownames(theta_ests2) <- field_names
    for(ff in field_names){
      for_f <- paste0(" for ",ff)
      hyper_f <- theta_ests[grep(for_f, names(theta_ests))]
      logtau_f <- theta_ests[paste0('log_tau', for_f)]
      logkappa_f <- theta_ests[paste0('log_kappa', for_f)]
      (range_f <- sqrt(8*nu)/exp(logkappa_f)) # r = sqrt(8*nu)/kappa
      (var_f <- gamma(nu)/(exp(logtau_f)^2 * (4*pi)^(d/2) * exp(logkappa_f)^(2*nu)))
      theta_ests2[ff,1] <- range_f
      theta_ests2[ff,2] <- var_f
    }

    RSS <- vector("list", nS)
    for (ss in seq(nS)) {
      RSS[[ss]] <- rowSums(t(matrix(
        as.matrix(c(BOLD[[ss]])) - (
          do.call(cbind, design[[ss]][valid_cols[ss,]]) %*% as.matrix(c(field_ests[[ss]]))
        ),
        nrow = nT[ss]
      ))^2)
    }
    # resids <- t(matrix(as.matrix(y) - (X %*% coefs), nrow = nT_ss))
    hyperpar_posteriors <- get_posterior_densities2(
      INLA_model_obj=INLA_model_obj, #spde,
      field_names
    ) #hyperparameter posterior densities

    #translate log_kappa to spatial range (cannot do the same for variance since it's a function of both kappa and tau)
    SPDEpar_posteriors1 <- hyperpar_posteriors[hyperpar_posteriors$param=='log_kappa',]
    SPDEpar_posteriors1$param <- 'SPDE_range'
    SPDEpar_posteriors1$value <- sqrt(8*nu)/exp(SPDEpar_posteriors1$value) #compute r = sqrt(8*nu)/kappa
    hyperpar_posteriors <- rbind(hyperpar_posteriors, SPDEpar_posteriors1)

    #construct object to be returned
    INLA_model_obj <- switch(return_INLA,
      trimmed=trim_INLA_model_obj(INLA_model_obj, minimal=FALSE),
      full=INLA_model_obj,
      minimal=trim_INLA_model_obj(INLA_model_obj, minimal=TRUE)
    )
    attr(INLA_model_obj, "format") <- return_INLA
  } else {
    field_ests <- lapply(result_classical, function(x){ x$estimates })
    attr(field_ests, "GLM_type") <- "classical"
  }

  # Tidy up and return. --------------------------------------------------------

  # Unmask.
  for (ss in seq(nS)) {
    field_ests[[ss]] <- unmask_Mdat2In(field_ests[[ss]], spatial$maskIn, spatial$maskMdat)

    if (!is.null(RSS)) {
      RSS[[ss]] <- unmask_Mdat2In(RSS[[ss]], spatial$maskIn, spatial$maskMdat)
    }

    result_classical[[ss]]$estimates <- unmask_Mdat2In(
      result_classical[[ss]]$estimates, spatial$maskIn, spatial$maskMdat)
    result_classical[[ss]]$SE_estimates <- unmask_Mdat2In(
      result_classical[[ss]]$SE_estimates, spatial$maskIn, spatial$maskMdat)
    result_classical[[ss]]$resids <- unmask_Mdat2In(
      result_classical[[ss]]$resids, spatial$maskIn, spatial$maskMdat)
    result_classical[[ss]]$RSS <- unmask_Mdat2In(
      result_classical[[ss]]$RSS, spatial$maskIn, spatial$maskMdat)
  }

  result <- list(
    field_estimates = field_ests, # new in 6.0: these are masked.
    INLA_model_obj = INLA_model_obj,
    hyperpar_posteriors = hyperpar_posteriors,
    theta_estimates = theta_ests, #kappa and tau
    theta_estimates2 = theta_ests2, #range and var
    RSS = RSS,
    result_classical = result_classical,
    #result_FIR = result_FIR,
    scrub_vec = scrub_vec,
    spatial = spatial,
    spde = spde,
    field_names = field_names,
    session_names = session_names,
    BOLD_QC = BOLD_QC,
    # For joint group model ~~~~~~~~~~~~~
    #posterior_Sig_inv = Sig_inv,
    y = y_all,
    X = XA_all_list,
    prewhiten_info = prewhiten_info,
    logkappa_vec = logkappa_vec,
    logtau_vec = logtau_vec,
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call = match.call()
  )
  class(result) <- "fit_bglm"

  result
}
