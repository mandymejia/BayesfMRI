#' Make Formula
#'
#' @param beta_names char vector of the names of each bbeta object in the environment
#' @param repl_names char vector of the names of each replicate object in the environment
#' @param model_name char name of the model for each beta (e.g. name of an spde object)
#' @param hyper_initial Optional vector of initial values for hyperparameters of each latent field OR a list with each element corresponding to one column of the X matrix
#'
#' @return A formula representing the Bayesian GLM to be passed to `inla()`
#' @export
#' @importFrom stats as.formula
#'
#' @examples \dontrun{}
make_formula <- function(beta_names, repl_names, model_name, hyper_initial=NULL){

	# Example:
	# beta_names = bbeta1, bbeta2, ...
	# repl_names = repl1, repl2, ...
	# model_name = spde
	# formula: y ~ -1 + f(bbeta1, model = spde, replicate = repl1) + f(bbeta2, model = spde_sh, replicate = repl2)

	# check length of beta_names, repl_names, hyper_initial

	nx <- length(beta_names)

	if(!is.null(hyper_initial)){
		#if hyper_list provided is a vector, repeat it nx times as a list
		if(!is.list(hyper_initial)){
			hyper_initial <- rep(list(hyper_initial), nx)
		}
		hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')
	} else {
		hyper_vec <- NULL
	}

	formula_vec <- paste0('f(',beta_names, ', model = ', model_name, ', replicate = ', repl_names, hyper_vec, ')')
	formula_vec <- c('y ~ -1', formula_vec)
	formula_str <- paste(formula_vec, collapse=' + ')
	return(as.formula(formula_str))
}
