#' Organize replicates
#'
#' beta and repl vectors are of length \eqn{nvox X n_sess X n_task}.
#' 	The ith repl vector is an indicator vector for the cells corresponding to the ith column of x.
#' 	The ith beta vector contains data indices (e.g. 1,...,V) in the cells corresponding to the ith column of x.
#'
#' @param n_sess The number of sessions sharing hyperparameters (can be different tasks)
#' @param beta_names Vector of names for each task
#' @inheritParams mesh_Param_either
#'
#' @return replicates vector and betas for sessions
#'
#' @keywords internal
#'
organize_replicates <- function(n_sess, beta_names, mesh){

  if(!(class(mesh) %in% c('inla.mesh','BayesfMRI.spde'))) stop('mesh must be of class inla.mesh  (for surface data, see `help(make_mesh)`) or BayesfMRI.spde (for subcortical data, see `help(create_spde_vol3D)`)')
  spatial <- mesh$idx$loc

  if(class(mesh) == "BayesfMRI.spde") { # This is for the subcortical case
    # We have additional locations within the mesh in this case, so we need to
    # account for them while also taking into account the fact that the
    # indices will increment between the regions
    # mesh$idx is a list of indices within each region
    spatial <- NULL # This will be a vector of indices across all regions within a group
    max_length <- 0 # Preset the number to add to each region
    for(r in 1:length(mesh$idx)) { # The length of mesh$idx gives the number of regions
      spatial <- c(spatial, mesh$idx[[r]] + max_length) # Add in the indices of region r
      max_length <- nrow(mesh$vertices[[r]]) + max_length # Gives the amount to increment
    }
  }

  nvox <- length(spatial)

	n_task <- length(beta_names)

	grps <- ((1:(n_sess*n_task) + (n_task-1)) %% n_task) + 1 # 1, 2, .. n_task, 1, 2, .. n_task, ...
	repls <- vector('list', n_task)
	betas <- vector('list', n_task)
	names(betas) <- beta_names
	for(i in 1:n_task){
		inds_i <- (grps == i)

		#set up replicates vectors
		sess_NA_i <- rep(NA, n_sess*n_task)
		sess_NA_i[inds_i] <- 1:n_sess
		repls[[i]] <- rep(sess_NA_i, each=nvox)
		names(repls)[i] <- paste0('repl',i)

		#set up ith beta vector with replicates for sessions
		NAs <- rep(NA, nvox)
		preNAs <- rep(NAs, times=(i-1))
		postNAs <- rep(NAs, times=(n_task-i))
		betas[[i]] <- rep(c(preNAs, spatial, postNAs), n_sess)
	}

	result <- list(betas=betas, repls=repls)
	return(result)

}
