#' Check INLA and PARDISO
#'
#' @param require_PARDISO Is PARDISO required? Default: \code{TRUE}.
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
check_INLA <- function(require_PARDISO=TRUE){

  # Check packages -------------------------------------------------------------

  # Check to see that the INLA package is installed
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("This function requires the `INLA` package. See www.r-inla.org/download")
  }

  # Check to see if PARDISO is installed
  if (require_PARDISO) {
    if (!exists("inla.pardiso.check", mode = "function")) {
      warning(paste(
        "Please update to the latest stable version of `INLA` for full functionality",
        "and `PARDISO` compatibility. See www.r-inla.org/download\n"
      ))
    } else {
      if (grepl("FAILURE", toupper(INLA::inla.pardiso.check()))) {
        warning("Consider enabling `PARDISO` for faster computation. See `inla.pardiso()`")
      } else {
        INLA::inla.setOption(smtp='pardiso')
      }
      #inla.pardiso()
    }
  }

  invisible(NULL)
}

#' Estimate INLA model
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param formula Formula to put into inla
#' @param data Dataset
#' @param A Large, sparse observation matrix
#' @param spde The spatial model, an object of class inla.spde
#' @param prec_initial Initial precision
#' @inheritParams num.threads_Param
#' @param int.strategy INLA strategy for numerical integration.  "eb" (empirical Bayes) is recommended for computational efficiency, or "ccd" for greater accuracy
#' @inheritParams verbose_Param_inla
#' @inheritParams contrasts_Param_inla
#' @param lincomb A linear combinations object created with \code{inla.make.lincomb}
#'
#' @return Results from INLA
#'
#' @export
estimate_model <- function(formula, data, A, spde, prec_initial, num.threads=4, int.strategy = "eb", verbose=FALSE, contrasts = NULL, lincomb = NULL){

  check_INLA(require_PARDISO=FALSE)

  INLA::inla(
    formula, data=data, control.predictor=list(A=A, compute = TRUE),
    verbose = verbose, keep = FALSE, num.threads = num.threads,
    control.inla = list(strategy = "gaussian", int.strategy = int.strategy),
    control.family=list(hyper=list(prec=list(initial=prec_initial))),
    control.compute=list(config=TRUE), contrasts = contrasts, lincomb = lincomb #required for excursions
  )
}


# Make Formula
#
# @param beta_names char vector of the names of each bbeta object in the environment
# @param repl_names char vector of the names of each replicate object in the environment
# @param hyper_initial Optional vector of initial values for hyperparameters of each latent field OR a list with each element corresponding to one column of the X matrix
#
# @return A formula representing the Bayesian GLM to be passed to `inla()`
#
# @importFrom stats as.formula
#
# @keywords internal
# make_formula <- function(beta_names, repl_names, hyper_initial=NULL){
#
#   # Example:
#   # beta_names = bbeta1, bbeta2, ...
#   # repl_names = repl1, repl2, ...
#   # formula: y ~ -1 + f(bbeta1, model = spde, replicate = repl1) + f(bbeta2, model = spde_sh, replicate = repl2)
#
#   # check length of beta_names, repl_names, hyper_initial
#
#   n_beta <- length(beta_names)
#
#   if(!is.null(hyper_initial)){
#     #if hyper_list provided is a vector, repeat it n_beta times as a list
#     if(!is.list(hyper_initial)){
#       hyper_initial <- rep(list(hyper_initial), n_beta)
#     }
#     hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')
#   } else {
#     hyper_vec <- NULL
#   }
#
#   formula_vec <- paste0('f(',beta_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
#   formula_vec <- c('y ~ -1', formula_vec)
#   formula_str <- paste(formula_vec, collapse=' + ')
#   return(formula_str)
# }

#' Make data list for \code{estimate_model}
#'
#' Make data list to be passed to \code{estimate_model}
#'
#' @param y Vectorized BOLD data (all voxels, sessions, etc.)
#' @param X List (length = number of sessions) of sparse design matrices size TVxVK from each session, each created using `organize_data()`
#' @param betas List (length = number of tasks) of bbeta objects from organize_replicates
#' @param repls List (length = number of tasks) of repl objects from organize_replicates
#'
#' @return List
#'
#' @importFrom Matrix bdiag
#'
#' @keywords internal
make_data_list <- function(y, X, betas, repls){

  # Check length/dimensions of y, X, elements of betas and repls all match
  n_sess <- length(X)
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

  return(model_data)
}


#' Extract Estimates of Activation
#'
#' Obtains the posterior mean or other summary statistic for each latent field
#'
#' @param object An object of class ‘"inla"’, a result of a call to inla
#' @param session_names Vector of fMRI session names
#' @param mask (Optional) Original mask applied to data before model fitting
#' @param stat A string representing the posterior summary statistic to be returned
#'
#' @return Estimates from inla model
#'
#' @keywords internal
extract_estimates <- function(object, session_names, mask=NULL, stat='mean'){

  if (!inherits(object, "inla")) { stop("Object is not of class 'inla'") }

  res.beta <- object$summary.random
  nbeta <- length(res.beta)
  beta_names <- names(res.beta)

  n_sess <- length(session_names)
  n_loc <- length(res.beta[[1]]$mean)/n_sess #number of locations for which beta is estimated
  if(!is.null(mask)) {
    V <- length(mask)
    if(sum(mask) != n_loc) warning('Number of nonzeros in mask does not equal the number of data locations in the model')
  } else {
    V <- n_loc
    mask <- rep(1,V)
  }
  betas <- vector('list', n_sess)
  names(betas) <- session_names

  stat_names <- names(res.beta[[1]])
  if(! (stat %in% stat_names) ) stop(paste0('stat must be one of following: ', paste(stat_names, collapse = ', ')))
  stat_ind <- which(stat_names==stat)


  for(v in 1:n_sess){
    inds_v <- (1:n_loc) + (v-1)*n_loc #indices of beta vector corresponding to session v
    betas_v <- matrix(NA, nrow=n_loc, ncol=nbeta)
    for(i in 1:nbeta){
      est_iv <- res.beta[[i]][[stat_ind]][inds_v]
      betas_v[,i] <- est_iv
    }
    betas[[v]] <- matrix(NA, nrow=V, ncol=nbeta)
    betas[[v]][mask==1,] <- betas_v
    colnames(betas[[v]]) <- beta_names
  }
  return(betas)
}


#' Extracts posterior density estimates for hyperparameters
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param object An object of class \code{"inla"}, a result of a call to
#'  \code{inla()}
#' @param spde The model used for the latent fields in the \code{inla()} call,
#'  an object of class \code{"inla.spde"}
#' @param beta_names Descriptive names of model regressors (tasks).
#'
#' @return Long-form data frame containing posterior densities for the
#'  hyperparameters associated with each latent field
#'
#' @keywords internal
get_posterior_densities <- function(object, spde, beta_names){

  numbeta <- length(beta_names)

  for(b in 1:numbeta){
    name_b <- beta_names[b]
    result.spde.b <- INLA::inla.spde2.result(object, name_b, spde)
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
