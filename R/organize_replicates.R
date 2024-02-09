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

	result <- list(betas=betas, repls=repls)
	return(result)

}
