#' Extract Estimates
#'
#' @param object An object of class ‘"inla"’, a result of a call to inla
#' @param mask Logical vector used to map beta estimates back to whole-brain field
#' @param session_names Name of session
#'
#' @return Estimates from inla model
#' @export
#'
#' @examples
extract_estimates <- function(object, mask, session_names){

  if(class(object) != "inla"){
    stop("Object is not of class 'inla'")
  }

	res.beta <- object$summary.random
	nbeta <- length(res.beta)
	beta_names <- names(res.beta)

	nvox2 <- sum(mask)
	n_sess <- length(session_names)

	if(length(res.beta[[1]]$mean) != nvox2*n_sess) stop('Length of estimate vectors must equal number of sessions x number of voxels')

	betas <- vector('list', n_sess)
	names(betas) <- session_names
	for(v in 1:n_sess){
		sess_name <- session_names[v]
		inds_v <- (1:nvox2) + (v-1)*nvox2 #indices of beta vector corresponding to session v
		betas_v <- matrix(rep(mask, nbeta), ncol=nbeta)

		for(i in 1:nbeta){
			est_iv <- res.beta[[i]]$mean[inds_v]
			betas_v[mask,i] <- est_iv
		}
		betas[[v]] <- betas_v
	}
	return(betas)
}

