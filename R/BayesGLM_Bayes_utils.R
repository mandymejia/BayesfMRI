#' Organize replicates
#'
#' beta and repl vectors are of length \eqn{n_mesh \times n_sess \times n_field}.
#' 	The ith repl vector is an indicator vector for the cells corresponding to the ith column of x.
#' 	The ith beta vector contains data indices (e.g. 1,...,V) in the cells corresponding to the ith column of x.
#'
#' @param n_sess The number of sessions sharing hyperparameters (can be different fields)
#' @param field_names Vector of names for each field
#' @param n_mesh Number of mesh locations
# @param data_loc Indices of original data locations
#'
#' @return replicates vector and betas for sessions
#'
#' @keywords internal
#'
organize_replicates <- function(n_sess, field_names, n_mesh){ #data_loc){

  spatial <- 1:n_mesh #data_loc
  #spatial <- mesh$idx$loc

	n_field <- length(field_names)

	grps <- ((1:(n_sess*n_field) + (n_field-1)) %% n_field) + 1 # 1, 2, .. n_field, 1, 2, .. n_field, ...
	repls <- vector('list', n_field)
	betas <- vector('list', n_field)
	names(betas) <- field_names
	for(i in 1:n_field){
		inds_i <- (grps == i)

		#set up replicates vectors
		sess_NA_i <- rep(NA, n_sess*n_field)
		sess_NA_i[inds_i] <- 1:n_sess
		repls[[i]] <- rep(sess_NA_i, each=n_mesh)
		names(repls)[i] <- paste0('repl',i)

		#set up ith beta vector with replicates for sessions
		NAs <- rep(NA, n_mesh)
		preNAs <- rep(NAs, times=(i-1))
		postNAs <- rep(NAs, times=(n_field-i))
		betas[[i]] <- rep(c(preNAs, spatial, postNAs), n_sess)
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
#' @param betas List (length = number of fields) of bbeta objects from organize_replicates
#' @param repls List (length = number of fields) of repl objects from organize_replicates
#'
#' @return List
#'
#' @importFrom Matrix bdiag
#'
#' @keywords internal
make_data_list <- function(y, X, betas, repls){

  # Check length/dimensions of y, X, elements of betas and repls all match
  nx <- length(betas) #check same as length(repls)
  #figure out nvox
  #check dim(X)
  #check length of betas and repls

  numel <- 1 + length(betas) + length(repls) + 1
  model_data <- vector('list', length=numel)
  names(model_data) <- c('y', 'X', names(betas), names(repls))
  model_data$y <- y
  model_data$X <- bdiag(X) #row/col structure: sess1_beta1, sess1_beta2, sess2_beta1, sess2_beta2, ...

  nbeta <- length(betas)
  for(i in 1:nbeta){
    model_data[[2+i]] <- betas[[i]]
    model_data[[2+nbeta+i]] <- repls[[i]]
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
#' @param mask (Optional) Original mask applied to data before model fitting
#' @param inds the indices in the mesh that correspond to the V data locations
#' @param stat A string representing the posterior summary statistic to be returned
#'
#' @return Estimates from inla model
#'
#' @keywords internal
extract_estimates <- function(INLA_model_obj, session_names, mask=NULL, inds, stat='mean'){

  if (!inherits(INLA_model_obj, "inla")) { stop("Object is not of class 'inla'") }

  res.beta <- INLA_model_obj$summary.random
  nbeta <- length(res.beta)
  field_names <- names(res.beta)

  nS <- length(session_names)

  #determine number of locations
  #	nV  = the number of data locations used in model fitting
  # nV2 or n_mesh = the number of mesh locations, which may be a superset
  # nV0 = the number of data locations prior to applying a mask
  if(is.null(mask)) mask <- rep(1, length(inds))
  nV0 <- length(mask) #number of data locations prior to applying a mask pre-model fitting
  nV <- length(inds) #number of data locations included in the model -- should match sum(mask)
  if(sum(mask) != nV) warning('Number of nonzeros in mask does not equal the number of data locations in the model')
  nV2 <- length(res.beta[[1]]$mean)/nS #number of locations in beta estimates

  betas <- vector('list', nS)
  names(betas) <- session_names

  stat_names <- names(res.beta[[1]])
  if(! (stat %in% stat_names) ) stop(paste0('stat must be one of following: ', paste(stat_names, collapse = ', ')))
  stat_ind <- which(stat_names==stat)

  for (ss in seq(nS)) {
    inds_ss <- (1:nV2) + (ss-1)*nV2 #indices of beta vector corresponding to session v
    betas_ss <- matrix(NA, nrow=nV2, ncol=nbeta)
    for (bb in seq(nbeta)) {
      est_iv <- res.beta[[bb]][[stat_ind]][inds_ss]
      betas_ss[,bb] <- est_iv
    }
    #remove boundary locations
    betas_ss <- betas_ss[inds,]
    #unmask, with NA for data locations excluded from analysis
    betas[[ss]] <- matrix(NA, nrow=nV0, ncol=nbeta)
    betas[[ss]][mask==1,] <- betas_ss
    colnames(betas[[ss]]) <- field_names
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
