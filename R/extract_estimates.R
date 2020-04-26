#' Extract Estimates of Activation
#'
#' @description Obtains the posterior mean or other summary statistic for each latent field
#'
#' @param object An object of class ‘"inla"’, a result of a call to inla
#' @param session_names Vector of fMRI session names
#' @param stat A string representing the posterior summary statistic to be returned
#'
#' @return Estimates from inla model
#' @export
#'
extract_estimates <- function(object, session_names, stat='mean'){

  if(class(object) != "inla"){
    stop("Object is not of class 'inla'")
  }

	res.beta <- object$summary.random
	nbeta <- length(res.beta)
	beta_names <- names(res.beta)

	n_sess <- length(session_names)
	n_loc <- length(res.beta[[1]]$mean)/n_sess
	betas <- vector('list', n_sess)
	names(betas) <- session_names

	stat_names <- names(res.beta[[1]])
	if(! (stat %in% stat_names) ) stop(paste0('stat must be one of following: ', paste(stat_names, collapse = ', ')))
	stat_ind <- which(stat_names==stat)


	for(v in 1:n_sess){
		inds_v <- (1:n_loc) + (v-1)*n_loc #indices of beta vector corresponding to session v
		betas_v <- matrix(NA, nrow=n_loc, ncol=nbeta)
		colnames(betas_v) <- beta_names

		for(i in 1:nbeta){
			est_iv <- res.beta[[i]][[stat_ind]][inds_v]
			betas_v[,i] <- est_iv
		}
		betas[[v]] <- betas_v
	}
	return(betas)
}

