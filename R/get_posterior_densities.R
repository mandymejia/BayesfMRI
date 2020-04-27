#' Extracts posterior density estimates for hyperparameters
#'
#' @param object An object of class ‘"inla"’, a result of a call to \code{inla()}
#' @param spde The model used for the latent fields in the \code{inla()} call, an object of class ‘"inla.spde"’
#' @param beta_names (Optional) Descriptive names of model regressors (tasks).
#'
#' @return Long-form data frame containing posterior densities for the hyperparameters associated with each latent field
#' @export
#' @importFrom INLA inla.spde2.result
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
get_posterior_densities <- function(object, spde, beta_names=NULL){

  beta_names_model <- names(object$summary.random)
  numbeta <- length(beta_names_model)
  if(numbeta != length(beta_names)) {
    warning('Length of beta_names invalid.  Setting to NULL.')
    beta_names <- NULL
  }

	for(b in 1:numbeta){
	  name_b = beta_names_model[b]
		result.spde.b = inla.spde2.result(object, name_b, spde)
		# Kappa and Tau
		log_kappa.b = as.data.frame(result.spde.b$marginals.log.kappa$kappa.1)
		log_tau.b = as.data.frame(result.spde.b$marginals.log.tau$tau.1)
		names(log_kappa.b) <- names(log_tau.b) <- c('value','density')
		log_kappa.b$param <- 'log_kappa'
		log_tau.b$param <- 'log_tau'
		df.b <- rbind(log_kappa.b, log_tau.b)
		df.b$beta <- name_b
		if(!is.null(beta_names)) df.b$name <- beta_names[b] else df.b$name <- NA
		if(b==beta_names_model[1]) df <- df.b else df <- rbind(df, df.b)
	}
	df <- df[,c('beta','name','param','value','density')]
	return(df)
}

