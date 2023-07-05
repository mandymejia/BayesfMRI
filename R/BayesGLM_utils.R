#' Check INLA and PARDISO
#'
#' @param require_PARDISO Is PARDISO required? Default: \code{FALSE}.
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
#' @param stat A string representing the posterior summary statistic to be returned
#'
#' @return Estimates from inla model
#'
#' @keywords internal
extract_estimates <- function(INLA_model_obj, session_names, mask=NULL, stat='mean'){

  if (!inherits(INLA_model_obj, "inla")) { stop("Object is not of class 'inla'") }

  res.beta <- INLA_model_obj$summary.random
  nbeta <- length(res.beta)
  task_names <- names(res.beta)

  nS <- length(session_names)
  n_loc <- length(res.beta[[1]]$mean)/nS #number of locations for which beta is estimated
  if(!is.null(mask)) {
    V <- length(mask)
    if(sum(mask) != n_loc) warning('Number of nonzeros in mask does not equal the number of data locations in the model')
  } else {
    V <- n_loc
    mask <- rep(1,V)
  }
  betas <- vector('list', nS)
  names(betas) <- session_names

  stat_names <- names(res.beta[[1]])
  if(! (stat %in% stat_names) ) stop(paste0('stat must be one of following: ', paste(stat_names, collapse = ', ')))
  stat_ind <- which(stat_names==stat)


  for (ss in seq(nS)) {
    inds_ss <- (1:n_loc) + (ss-1)*n_loc #indices of beta vector corresponding to session v
    betas_ss <- matrix(NA, nrow=n_loc, ncol=nbeta)
    for (bb in seq(nbeta)) {
      est_iv <- res.beta[[bb]][[stat_ind]][inds_ss]
      betas_ss[,bb] <- est_iv
    }
    betas[[ss]] <- matrix(NA, nrow=V, ncol=nbeta)
    betas[[ss]][mask==1,] <- betas_ss
    colnames(betas[[ss]]) <- task_names
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
#' @param task_names Descriptive names of model regressors (tasks).
#'
#' @return Long-form data frame containing posterior densities for the
#'  hyperparameters associated with each latent field
#'
#' @keywords internal
get_posterior_densities <- function(INLA_model_obj, spde, task_names){

  numbeta <- length(task_names)

  for(b in 1:numbeta){
    name_b <- task_names[b]
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
    tasks = object$task_names,
    sessions = object$session_names,
    n_sess_orig = object$n_sess_orig,
    n_loc_total = length(object$mask),
    n_loc_modeled = sum(object$mask),
    GLM_type = attr(object$task_estimates, "GLM_type")
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
  cat("Tasks:    ", paste0("(", length(x$tasks), ") ", paste(x$tasks, collapse=", ")), "\n")
  if (length(x$sessions)==1 && x$sessions == "session_combined") {
    cat("Sessions: ", paste0("(", x$n_sess_orig, ", combined) \n"))
  } else {
    cat("Sessions: ", paste0("(", length(x$sessions), ") ", paste(x$sessions, collapse=", ")), "\n")
  }
  cat("Locations:", x$n_loc_modeled, "modeled,", x$n_loc_total, "total", "\n")
  cat("GLM type: ", x$GLM_type, "\n")
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
    tasks = x[[1]]$tasks,
    sessions = x[[1]]$sessions,
    n_sess_orig = x[[1]]$n_sess_orig,
    n_loc_total = lapply(x, '[[', "n_loc_total"),
    n_loc_modeled = lapply(x, '[[', "n_loc_modeled"),
    #xii = summary(x$task_estimates_xii$classical[[1]]),
    GLM_type = x[[1]]$GLM_type
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
  cat("Tasks:    ", paste0("(", length(x$tasks), ") ", paste(x$tasks, collapse=", ")), "\n")
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
  cat("GLM type: ", x$GLM_type, "\n")
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
