#' Make replicates
#'
#' beta and repl vectors are of length \eqn{nMesh \times nSess \times n_field}.
#' 	The ith repl vector is an indicator vector for the cells corresponding to the ith column of x.
#' 	The ith beta vector contains data indices (e.g. 1,...,V) in the cells corresponding to the ith column of x.
#'
#' @param nSess The number of sessions sharing hyperparameters (can be different fields)
#' @param field_names Vector of names for each field
#' @param nMesh Number of mesh locations
# @param data_loc Indices of original data locations
#'
#' @return replicates vector and betas for sessions
#'
#' @keywords internal
#'
make_replicates <- function(nSess, field_names, nMesh){ #data_loc){

  seq_nMesh <- seq(nMesh) #data_loc
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
		repls[[ii]] <- rep(sess_NA_ii, each=nMesh)
		names(repls)[ii] <- paste0('repl',ii)

		#set up ith beta vector with replicates for sessions
		NAs <- rep(NA, nMesh)
		preNAs <- rep(NAs, times=(ii-1))
		postNAs <- rep(NAs, times=(nK-ii))
		betas[[ii]] <- rep(c(preNAs, seq_nMesh, postNAs), nSess)
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
#' @param X List (length = number of sessions) of sparse design matrices size TVxVK from each session, each created using `sparse_and_PW()`
#' @param betas List (length = number of fields) of bbeta objects from make_replicates
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
  spatial, spatial_type, spde, stat='mean'){

  if (!inherits(INLA_model_obj, "inla")) { stop("Object is not of class 'inla'") }

  res.beta <- INLA_model_obj$summary.random
  nK <- length(res.beta)
  field_names <- names(res.beta)

  nS <- length(session_names)
  mask <- if (spatial_type=="mesh") {
    spatial$mask
  } else {
    spatial$labels != 0
  }

  #determine number of locations
  #	nV  = the number of data locations used in model fitting
  # nV2 or nMesh = the number of mesh locations, which may be a superset
  # nV0 = the number of data locations prior to applying a mask
  nV_T <- length(mask) #number of data locations prior to applying a mask pre-model fitting
  nV_DB <- if (spatial_type=="mesh") {
    sum(mask)
  } else {
    get_nV(spatial, "voxel")$DB
  }
  stopifnot(nV_DB*nS == length(res.beta[[1]]$mean))

  betas <- vector('list', nS)
  names(betas) <- session_names

  stat_names <- names(res.beta[[1]])
  if(! (stat %in% stat_names) ) stop(paste0('stat must be one of following: ', paste(stat_names, collapse = ', ')))
  stat_ind <- which(stat_names==stat)

  for (ss in seq(nS)) {
    inds_ss <- seq(nV_DB) + (ss-1)*nV_T
    betas[[ss]] <- do.call(cbind,
      lapply(setNames(seq(nK), field_names), function(kk){
        res.beta[[kk]][[stat_ind]][inds_ss]
      })
    )
    if (spatial_type=="voxel") {
      betas[[ss]] <- betas[[ss]][spatial$buffer_mask,,drop=FALSE]
    }
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
    which_kappa_b <- which_b[2] #Theta2
    which_tau_b <- which_b[1] #Theta1
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
#'  \code{id_activations} or \code{BayesGLM2}.
#'
#' @param INLA_model_obj An object of class \code{"inla"}.
#' @param minimal Just keep the two parameters needed for \code{BayesGLM2}?
#'  Default: \code{FALSE}. \code{!minimal} is required for
#'  \code{id_activations}, but \code{minimal} is sufficient for
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
