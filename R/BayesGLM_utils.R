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
#' @param X List (length = number of sessions) of sparse design matrices size TVxVK from each session, each created using `organize_data()`
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

  return(model_data)
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
  return(betas)
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

#' Summarize a \code{"BayesGLM"} object
#'
#' Summary method for class \code{"BayesGLM"}
#'
#' @param object Object of class \code{"BayesGLM"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BayesGLM"} object, a list summarizing the properties
#'  of \code{object}.
#' @method summary BayesGLM
summary.BayesGLM <- function(object, ...) {

  x <- list(
    fields = object$field_names,
    sessions = object$session_names,
    n_sess_orig = object$n_sess_orig,
    n_loc_total = length(object$mask),
    n_loc_modeled = sum(object$mask),
    GLM_type = attr(object$field_estimates, "GLM_type"),
    design_is_multiple = !is.null(object$result_multiple)
  )

  class(x) <- "summary.BayesGLM"

  return(x)
}

#' @rdname summary.BayesGLM
#' @export
#'
#' @param x Object of class \code{"summary.BayesGLM"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BayesGLM
print.summary.BayesGLM <- function(x, ...) {
  cat("====BayesGLM result====================\n")
  cat("Fields:   ", paste0("(", length(x$fields), ") ", paste(x$fields, collapse=", ")), "\n")
  if (length(x$sessions)==1 && x$sessions == "session_combined") {
    cat("Sessions: ", paste0("(", x$n_sess_orig, ", combined) \n"))
  } else {
    cat("Sessions: ", paste0("(", length(x$sessions), ") ", paste(x$sessions, collapse=", ")), "\n")
  }
  cat("Locations:", x$n_loc_modeled, "modeled,", x$n_loc_total, "total", "\n")
  if (x$design_is_multiple) {
    cat("GLM type: ", "classical (multiple designs)", "\n")
  } else {
    cat("GLM type: ", x$GLM_type, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.BayesGLM
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BayesGLM
print.BayesGLM <- function(x, ...) {
  print.summary.BayesGLM(summary(x))
}

#' Summarize a \code{"BayesGLM_cifti"} object
#'
#' Summary method for class \code{"BayesGLM_cifti"}
#'
#' @param object Object of class \code{"BayesGLM_cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BayesGLM_cifti"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary BayesGLM_cifti
summary.BayesGLM_cifti <- function(object, ...) {

  x <- lapply(object$BayesGLM_results, summary)
  x <- x[!vapply(object$BayesGLM_results, is.null, FALSE)]
  x <- list(
    fields = x[[1]]$fields,
    sessions = x[[1]]$sessions,
    n_sess_orig = x[[1]]$n_sess_orig,
    n_loc_total = lapply(x, '[[', "n_loc_total"),
    n_loc_modeled = lapply(x, '[[', "n_loc_modeled"),
    #xii = summary(x$estimates_xii$classical[[1]]),
    GLM_type = x[[1]]$GLM_type,
    design_is_multiple = x[[1]]$design_is_multiple
  )
  class(x) <- "summary.BayesGLM_cifti"

  return(x)
}

#' @rdname summary.BayesGLM_cifti
#' @export
#'
#' @param x Object of class \code{"summary.BayesGLM_cifti"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BayesGLM_cifti
print.summary.BayesGLM_cifti <- function(x, ...) {
  cat("====BayesGLM_cifti result==============\n")
  cat("Fields:   ", paste0("(", length(x$fields), ") ", paste(x$fields, collapse=", ")), "\n")
  if (length(x$sessions)==1 && x$sessions == "session_combined") {
    cat("Sessions: ", paste0("(", x$n_sess_orig, ", combined) \n"))
  } else {
    cat("Sessions: ", paste0("(", length(x$sessions), ") ", paste(x$sessions, collapse=", ")), "\n")
  }
  cat("Locations:\n")
  for (ii in seq(length(x$n_loc_total))) {
    cat(
      "          ", paste0(names(x$n_loc_total)[ii], ": ", x$n_loc_modeled[[ii]]),
      "modeled,", x$n_loc_total[[ii]], "total", "\n"
    )
  }
  if (x$design_is_multiple) {
    cat("GLM type: ", "classical (multiple designs)", "\n")
  } else {
    cat("GLM type: ", x$GLM_type, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.BayesGLM_cifti
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BayesGLM_cifti
print.BayesGLM_cifti <- function(x, ...) {
  print.summary.BayesGLM_cifti(summary(x))
}
