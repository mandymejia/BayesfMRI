#' Extracts posterior density estimates for hyperparameters
#'
#' @param object An object of class ‘"inla"’, a result of a call to \code{inla()}
#' @param spde The model used for the latent fields in the \code{inla()} call, an object of class ‘"inla.spde"’
#'
#' @return Long-form data frame containing posterior densities for the hyperparameters associated with each latent field
#' @export
#' @importFrom INLA inla.spde2.result
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
get_posterior_densities <- function(object, spde){

  beta_names <- names(object$summary.random)

	for(b in beta_names){
		result.spde.b = inla.spde2.result(object, b, spde)
		# Kappa and Tau
		log_kappa.b = as.data.frame(result.spde.b$marginals.log.kappa$kappa.1)
		log_tau.b = as.data.frame(result.spde.b$marginals.log.tau$tau.1)
		names(log_kappa.b) <- names(log_tau.b) <- c('value','density')
		log_kappa.b$param <- 'log_kappa'
		log_tau.b$param <- 'log_tau'
		df.b <- rbind(log_kappa.b, log_tau.b)
		df.b$beta <- b
		if(b==beta_names[1]) df <- df.b else df <- rbind(df, df.b)
	}
	df <- df[,c('beta','param','value','density')]
	return(df)
}

