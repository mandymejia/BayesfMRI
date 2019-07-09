#' Organize replicates
#'
#' @param n_sess The number of sessions sharing hyperparameters (can be different tasks)
#' @param n_task Number of regressors or tasks
#' @param mesh The mesh for the data
#'
#' @return replicates vector and betas for sessions
#' @export
#'
#' @examples \dontrun{}
organize_replicates <- function(n_sess, n_task, mesh){

	# create vectors for each task:
	#

	# beta and repl vectors are of length nvox * n_sess * n_task
	# ith repl vector is an indicator vector for the cells corresponding to the ith column of x
	# ith beta vector contains data indices (e.g. 1,...,V) in the cells corresponding to the ith column of x

	spatial <- mesh$idx$loc
	nvox <- length(spatial)

	grps <- ((1:(n_sess*n_task) + (n_task-1)) %% n_task) + 1 # 1, 2, .. n_task, 1, 2, .. n_task, ...
	repls <- vector('list', n_task)
	betas <- vector('list', n_task)
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
		names(betas)[i] <- paste0('bbeta',i)
	}

	result <- list(betas=betas, repls=repls)
	return(result)

}
