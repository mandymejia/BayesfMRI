#' Make replicates
#'
#' beta and repl vectors are of length \eqn{nMesh \times nSess \times n_field}.
#' 	The ith repl vector is an indicator vector for the cells corresponding to the ith column of x.
#' 	The ith beta vector contains data indices (e.g. 1,...,V) in the cells corresponding to the ith column of x.
#'
#' @param nSess The number of sessions sharing hyperparameters (can be different fields)
#' @param field_names Vector of names for each field
#' @param spatial Spatial info
#'
#' @return replicates vector and betas for sessions
#'
#' @keywords internal
#'
make_replicates <- function(nSess, field_names, spatial){

  nV <- get_nV(spatial)
	nK <- length(field_names)

	grps <- ((1:(nSess*nK) + (nK-1)) %% nK) + 1 # 1, 2, .. nK, 1, 2, .. nK, ...

  betas <- repls <- vector('list', nK)
  names(betas) <- field_names
  names(repls) <- paste0("repl", seq(nK))

	for (ii in 1:nK) {
		inds_ii <- (grps == ii)

		#set up replicates vectors
		sess_NA_ii <- rep(NA, nSess*nK)
		sess_NA_ii[inds_ii] <- 1:nSess
		repls[[ii]] <- rep(sess_NA_ii, each=nV$model)
		names(repls)[ii] <- paste0('repl',ii)

		#set up ith beta vector with replicates for sessions
		NAs <- rep(NA, nV$model)
		preNAs <- rep(NAs, times=(ii-1))
		postNAs <- rep(NAs, times=(nK-ii))
		betas[[ii]] <- rep(c(preNAs, seq(nV$model), postNAs), nSess)
	}

	list(betas=betas, repls=repls)
}

#' Check INLA and PARDISO
#'
#' @param require_PARDISO Is PARDISO required? Default: \code{FALSE}.
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
check_INLA <- function(require_PARDISO=FALSE){

  # Check packages -------------------------------------------------------------

  # Check to see that the INLA package is installed
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("This function requires the `INLA` package. See www.r-inla.org/download")
  }

  # Check to see if PARDISO is installed
  if (require_PARDISO) {
    if (any(grepl("FAILURE", toupper(INLA::inla.pardiso.check())))) {
      warning("Consider enabling `PARDISO` for faster computation. See `inla.pardiso()`")
    } else {
      INLA::inla.setOption(smtp='pardiso')
    }
    #inla.pardiso()
  }

  invisible(NULL)
}

#' Make data list for \code{estimate_model}
#'
#' Make data list to be passed to \code{estimate_model}
#'
#' @param y Vectorized BOLD data (all voxels, sessions, etc.)
#' @param X List (length = number of sessions) of sparse design matrices size TVxVK from each session, each created using \code{sparse_and_PW}
#' @param betas List (length = number of fields) of beta objects from make_replicates
#' @param repls List (length = number of fields) of repl objects from make_replicates
#'
#' @return List
#'
#' @importFrom Matrix bdiag
#'
#' @keywords internal
make_data_list <- function(y, X, betas, repls){

  # Check length/dimensions of y, X, elements of betas and repls all match
  nK <- length(betas) #check same as length(repls)
  #figure out nvox
  #check dim(X)
  #check length of betas and repls

  numel <- 1 + length(betas) + length(repls) + 1
  model_data <- vector('list', length=numel)
  names(model_data) <- c('y', 'X', names(betas), names(repls))
  model_data$y <- y
  model_data$X <- bdiag(X) #row/col structure: sess1_beta1, sess1_beta2, sess2_beta1, sess2_beta2, ...

  for (kk in seq(nK)) {
    model_data[[2+kk]] <- betas[[kk]]
    model_data[[2+nK+kk]] <- repls[[kk]]
  }

  model_data
}

#' Extract Estimates of Activation
#'
#' Obtains the posterior mean or other summary statistic for each latent field
#'
#' @param INLA_model_obj An object of class \code{"inla"}, a result of a call to
#'  \code{inla}.
#' @param session_names Vector of fMRI session names
#' @param stat A string representing the posterior summary statistic to be returned
#'
#' @return Estimates from inla model
#'
#' @keywords internal
extract_estimates <- function(
  INLA_model_obj, session_names,
  spatial, spde, stat='mean'){

  if (!inherits(INLA_model_obj, "inla")) { stop("Object is not of class 'inla'") }

  res.beta <- INLA_model_obj$summary.random
  nK <- length(res.beta)
  field_names <- names(res.beta)
  nS <- length(session_names)

  nV <- get_nV(spatial)

  stopifnot(nV$model*nS == length(res.beta[[1]]$mean))

  betas <- vector('list', nS)
  names(betas) <- session_names

  stat_names <- names(res.beta[[1]])
  if(! (stat %in% stat_names) ) stop(paste0('stat must be one of following: ', paste(stat_names, collapse = ', ')))
  stat_ind <- which(stat_names==stat)

  for (ss in seq(nS)) {
    inds_ss <- seq(nV$model) + (ss-1)*nV$total
    betas[[ss]] <- do.call(cbind,
      lapply(setNames(seq(nK), field_names), function(kk){
        res.beta[[kk]][[stat_ind]][inds_ss]
      })
    )
    betas[[ss]] <- betas[[ss]][spatial$Mmap,,drop=FALSE]
  }

  attr(betas, "GLM_type") <- "Bayesian"
  betas
}


#' Extracts posterior density estimates for hyperparameters
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param INLA_model_obj An object of class \code{"inla"}, a result of a call to
#'  \code{inla()}
#' @param spde The model used for the latent fields in the \code{inla()} call,
#'  an object of class \code{"inla.spde"}
#' @param field_names Descriptive names of model regressors (fields).
#'
#' @return Long-form data frame containing posterior densities for the
#'  hyperparameters associated with each latent field
#'
#' @keywords internal
get_posterior_densities <- function(INLA_model_obj, spde, field_names){

  numbeta <- length(field_names)

  for(b in 1:numbeta){
    name_b <- field_names[b]
    result.spde.b <- INLA::inla.spde2.result(INLA_model_obj, name_b, spde)
    # Kappa and Tau
    log_kappa.b <- as.data.frame(result.spde.b$marginals.log.kappa$kappa.1)
    log_tau.b <- as.data.frame(result.spde.b$marginals.log.tau$tau.1)
    names(log_kappa.b) <- names(log_tau.b) <- c('value','density')
    log_kappa.b$param <- 'log_kappa'
    log_tau.b$param <- 'log_tau'
    df.b <- rbind(log_kappa.b, log_tau.b)
    df.b$beta <- name_b
    if(b == 1) df <- df.b else df <- rbind(df, df.b)
  }

  df[,c('beta','param','value','density')]
}

#' Extracts posterior density estimates for hyperparameters for volumetric SPDE
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param INLA_model_obj An object of class \code{"inla"}, a result of a call to
#'  \code{inla()}
#' @param field_names Descriptive names of model regressors (fields).
#'
#' @return Long-form data frame containing posterior densities for the
#'  hyperparameters associated with each latent field
#'
#' @keywords internal
get_posterior_densities2 <- function(INLA_model_obj, field_names){

  #marginal densities of hyperparameters
  hyper <- INLA_model_obj$marginals.hyperpar
  names_hyper <- names(INLA_model_obj$marginals.hyperpar)

  #grab marginals for kappa and tau for each latent field
  numbeta <- length(field_names)
  for(b in 1:numbeta){
    name_b <- field_names[b]
    which_b <- grep(name_b, names_hyper)
    which_kappa_b <- which_b[2] #Theta2 = log(kappa)
    which_tau_b <- which_b[1] #Theta1 = log(tau)
    log_kappa.b <- as.data.frame(hyper[[which_kappa_b]])
    log_tau.b <- as.data.frame(hyper[[which_tau_b]])
    names(log_kappa.b) <- names(log_tau.b) <- c('value','density')
    log_kappa.b$param <- 'log_kappa'
    log_tau.b$param <- 'log_tau'
    df.b <- rbind(log_kappa.b, log_tau.b)
    df.b$beta <- name_b
    if(b == 1) df <- df.b else df <- rbind(df, df.b)
  }

  #grab marginal for the precision of the likelihood
  df.prec <- as.data.frame(hyper[[1]])
  names(df.prec) <- c('value','density')
  df.prec$param <- 'prec'
  df.prec$beta <- 'none'
  df <- rbind(df, df.prec)

  df[,c('beta','param','value','density')]
}

#' Trim INLA object
#'
#' Trim an INLA object to only include what is necessary for
#'  \code{activations} or \code{BayesGLM2}.
#'
#' @param INLA_model_obj An object of class \code{"inla"}.
#' @param minimal Just keep the two parameters needed for \code{BayesGLM2}?
#'  Default: \code{FALSE}. \code{!minimal} is required for
#'  \code{activations}, but \code{minimal} is sufficient for
#'  \code{BayesGLM2}.
#'
#' @return A trimmed \code{"inla"} object.
#' @keywords internal
trim_INLA_model_obj <- function(INLA_model_obj, minimal=FALSE) {
  if (!inherits(INLA_model_obj, "inla")) {
    stop("This function only applies to objects with the 'inla' class.")
  }

  out_object <- list()
  out_object$misc$theta.mode <- INLA_model_obj$misc$theta.mode
  out_object$misc$cov.intern <- INLA_model_obj$misc$cov.intern
  if (!minimal) {
    out_object$.args$control.compute$config <- INLA_model_obj$.args$control.compute$config
    out_object$marginals.random <- INLA_model_obj$marginals.random
    out_object$misc$configs <- INLA_model_obj$misc$configs
  }

  class(out_object) <- "inla"
  out_object
}

#' Make \code{log_kappa} and \code{log_tau}
#'
#' Make \code{log_kappa} and \code{log_tau}
#'
#' @param spatial,hyperpriors,verbose See \code{fit_bayesglm}
#'
#' @keywords internal
log_kappa_tau <- function(spatial, hyperpriors, verbose){

  d <- switch(spatial$spatial_type, vertex=2, voxel=3)
  nu <- 2 - d/2 #nu = alpha - d/2, alpha = 2

  # Log kappa ----------------------------------------------------------------
  #determine reasonable values for spatial range
  if (spatial$spatial_type == "vertex") {
    ssv <- spatial$surf$vertices
    ssf <- spatial$surf$faces
    max_dist <- apply(ssv, 2, function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) #MAX distance within the mesh (Euclidean, a bad proxy for geodesic)
    range2 <- max(max_dist)/2 #this is our upper limit for the spatial correlation range
    tri_dist <- apply(
      matrix(1:nrow(ssf), ncol=1), 1, #for each face/triangle #MIN distance within the mesh
      function(x) { mean(dist(ssv[ssf[x,],])) } #compute the avg distance among vertices in the triangle
    )
    min_dist <- min(tri_dist) #distance between vertices in the smallest triangle
    range1 <- min_dist*2
    range0 <- range1*5

  } else {
    res <- abs(diag(spatial$trans_mat)[1:3])
    range2 <- c()
    for(r in levels(spatial$labels)){
      #get mask of just this ROI
      mask_r <- spatial$maskIn
      mask_r[mask_r][spatial$labels != r] <- FALSE
      #compute max distance within mask in each direction
      x_r <- diff(range(which(apply(mask_r, 1, sum) > 0)))*res[1]
      y_r <- diff(range(which(apply(mask_r, 2, sum) > 0)))*res[2]
      z_r <- diff(range(which(apply(mask_r, 3, sum) > 0)))*res[3]
      range2 <- c(range2, max(c(x_r, y_r, z_r))) #largest 1D distance in any direction
    }
    #this is our upper limit for the spatial correlation range. Since volumetric structures are smaller, we allow a larger range
    range2 <- max(range2)*2 #max over ROIs, doubled (in case we don't observe the full range)
    range1 <- min(res)*2 #smallest voxel dimension
    #reasonable starting value
    range0 <- range1*5
  }

  range_vec <- c(range1, range2, range0)

  if(hyperpriors == "informative") logkappa_vec <- log(sqrt(8*nu)/range_vec) #r = sqrt(8*nu)/kappa
  logkappa0 <- log(sqrt(8*nu)/range0) #we will use a good starting value even with the default prior

  # logkappa1 <-
  # logkappa2 <- log(sqrt(8*nu)/range2)
  # logkappa0 <- log(sqrt(8*nu)/range0) #starting value


  # Log tau ------------------------------------------------------------------
  # Get initial value for tau
  var2logtau <- function(var, d, kappa){
    nu <- 2 - d/2 #nu = alpha - d/2, alpha = 2
    tausq <- gamma(nu)/(var * (4*pi)^(d/2) * kappa^(2*nu))
    log(sqrt(tausq))
  }

  # # Get initial value for tau based on variance of classical GLM
  # var0 <- apply(result_classical[[ss]]$estimates, 2, var, na.rm=TRUE)

  var0 <- 0.1 #we usually expect most of the activation amplitudes to be between (-2 and 2) --> SD ~= 0.33, Var ~= 0.1
  var_vec <- c(0.01, 1, var0) #reasonable range for variance
  if(hyperpriors == "informative") logtau_vec <- var2logtau(var_vec, d, exp(logkappa0))
  logtau0 <- var2logtau(var0, d, exp(logkappa0)) #we will use a good starting value even with the default prior

  ### SUMMARY

  if (verbose > 0 && hyperpriors == 'informative') {
    cat(paste0('\tPutting an informative prior on kappa so that the spatial range is between ',
      round(range1, 2), ' and ', round(range2, 2), ' mm.\n',
      '\t\tLog kappa prior range (95% density): ', round(logkappa_vec[2], 2), ' to ', round(logkappa_vec[1], 2), '\n'))
    cat(paste0('\tPutting an informative prior on tau so variance of the spatial field is between ',
    (var_vec[1]), ' and ', (var_vec[2]), '.\n',
    '\t\tLog tau prior range (95% density): ',round(logtau_vec[2], 2), ' to ', round(logtau_vec[1], 2), '\n'))
  }

  list(
    logkappa_vec=logkappa_vec,
    logkappa0=logkappa0,
    logtau0=logtau0,
    logtau_vec=logtau_vec
  )
}
