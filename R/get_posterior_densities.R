#' Extracts posterior density estimates for hyperparameters
#'
#' @param result Result
#' @param spde SPDE
#' @param beta_names Names of covariates
#'
#' @return Long-form data frame
#' @export
#'
#' @examples \dontrun{}
get_posterior_densities <- function(result, spde, beta_names){

	for(b in beta_names){
		result.spde.b = inla.spde2.result(result, b, spde)

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

