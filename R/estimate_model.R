#' Estimate INLA model
#'
#' @param formula Formula to put into inla
#' @param data Dataset
#' @param A Large, sparse observation matrix
#' @param prec_initial Initial precision
#' @param num.threads Number of threads
#' @param int.strategy INLA strategy for numerical integration.  "eb" (empirical Bayes) is recommended for computational efficiency, or "ccd" for greater accuracy
#'
#' @return Results from INLA
#' @export
#'
#' @examples
estimate_model <- function(formula, data, A, prec_initial, num.threads=4, int.strategy = "eb"){

	require(INLA)
	result <- inla(formula, data=data, control.predictor=list(A=A, compute = TRUE),
					verbose = TRUE, keep = FALSE, num.threads = num.threads,
					control.inla = list(strategy = "gaussian", int.strategy = int.strategy),
					control.family=list(hyper=list(prec=list(initial=prec_initial))),
	    			control.compute=list(config=TRUE))
	return(result)
}
