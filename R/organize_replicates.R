#' Organize replicates
#'
#' beta and repl vectors are of length \eqn{nvox \times n_sess \times n_task}.
#' 	The ith repl vector is an indicator vector for the cells corresponding to the ith column of x.
#' 	The ith beta vector contains data indices (e.g. 1,...,V) in the cells corresponding to the ith column of x.
#'
#' @param n_sess The number of sessions sharing hyperparameters (can be different tasks)
#' @param task_names Vector of names for each task
#' @inheritParams mesh_Param_inla
# @inheritParams mesh_Param_either
#'
#' @return replicates vector and betas for sessions
#'
#' @keywords internal
#'
organize_replicates <- function(n_sess, task_names, mesh){

  if (!(inherits(mesh, "inla.mesh") || inherits(mesh, "BayesfMRI.spde"))) {
	stop('mesh must be of class inla.mesh  (for surface data, see `help(make_mesh)`) or BayesfMRI.spde (for subcortical data, see `help(create_spde_vol3D)`)')
  }
  spatial <- mesh$idx$loc

	nvox <- length(spatial)

	n_task <- length(task_names)

	grps <- ((1:(n_sess*n_task) + (n_task-1)) %% n_task) + 1 # 1, 2, .. n_task, 1, 2, .. n_task, ...
	repls <- vector('list', n_task)
	betas <- vector('list', n_task)
	names(betas) <- task_names
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
