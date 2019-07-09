#' Make data list
#'
#' @param y Vectorized BOLD data (all voxels, sessions, etc.)
#' @param X List of sparse design matrices from each session, each created using organize_data
#' @param betas List of bbeta objects from organize_replicates
#' @param repls List of repl objects from organize_replicates
#'
#' @return List
#'
#' @examples \dontrun{}
make_data_list <- function(y, X, betas, repls){

	# Check length/dimensions of y, X, elements of betas and repls all match
	n_sess <- length(X_all_list)
	nx <- length(betas) #check same as length(repls)
	#figure out nvox
	#check dim(X)
	#check length of betas and repls

	numel <- 1 + length(betas) + length(repls) + 1
	model_data <- vector('list', length=numel)
	names(model_data) <- c('y', 'X', names(betas), names(repls))
	model_data$y <- y
	model_data$X <- bdiag(X) #row/col structure: task1_sess1_beta1, task1_sess1_beta2, task1_sess2_beta1, task1_sess2_beta2, ..., task2_sess1_beta1, ...

	nbeta <- length(betas)
	for(i in 1:nbeta){
		model_data[[2+i]] <- betas[[i]]
		model_data[[2+nbeta+i]] <- repls[[i]]
	}

	return(model_data)
}
