#' Make Formula
#'
#' @param beta_names char vector of the names of each bbeta object in the environment
#' @param repl_names char vector of the names of each replicate object in the environment
#' @param spde The spatial model for the latent fields
#' @param hyper_initial Optional vector of initial values for hyperparameters of each latent field OR a list with each element corresponding to one column of the X matrix
#'
#' @return A formula representing the Bayesian GLM to be passed to `inla()`
#' @importFrom stats as.formula
#'
#' @examples \dontrun{}
make_formula <- function(beta_names, repl_names, spde, hyper_initial=NULL){

	# Example:
	# beta_names = bbeta1, bbeta2, ...
	# repl_names = repl1, repl2, ...
	# formula: y ~ -1 + f(bbeta1, model = spde, replicate = repl1) + f(bbeta2, model = spde_sh, replicate = repl2)

	# check length of beta_names, repl_names, hyper_initial

	n_beta <- length(beta_names)

	if(!is.null(hyper_initial)){
		#if hyper_list provided is a vector, repeat it n_beta times as a list
		if(!is.list(hyper_initial)){
			hyper_initial <- rep(list(hyper_initial), n_beta)
		}
		hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')
	} else {
		hyper_vec <- NULL
	}

	formula_vec <- paste0('f(',beta_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
	formula_vec <- c('y ~ -1', formula_vec)
	formula_str <- paste(formula_vec, collapse=' + ')
	formula <- as.formula(formula_str)
	environment(formula) <- globalenv()
	return(formula)
}
