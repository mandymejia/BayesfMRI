#' Estimate INLA model
#'
#' @param formula Formula to put into inla
#' @param data Dataset
#' @param A Large, sparse observation matrix
#' @param spde The spatial model, an object of class inla.spde
#' @param prec_initial Initial precision
#' @param num.threads Number of threads
#' @param int.strategy INLA strategy for numerical integration.  "eb" (empirical Bayes) is recommended for computational efficiency, or "ccd" for greater accuracy
#' @param verbose Logical indicating if should run in a verbose mode (default FALSE).
#'
#' @return Results from INLA
#' @export
#' @importFrom INLA inla
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
estimate_model <- function(formula, data, A, spde, prec_initial, num.threads=4, int.strategy = "eb", verbose=FALSE){

	result <- inla(formula, data=data, control.predictor=list(A=A, compute = TRUE),
					verbose = verbose, keep = FALSE, num.threads = num.threads,
					control.inla = list(strategy = "gaussian", int.strategy = int.strategy),
					control.family=list(hyper=list(prec=list(initial=prec_initial))),
	    			control.compute=list(config=TRUE)) #required for excursions
	return(result)
}
