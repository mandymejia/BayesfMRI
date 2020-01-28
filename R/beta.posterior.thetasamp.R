#' Internal function used in joint approach to group-analysis
#'
#' @param theta A sample of theta (hyperparameters)
#' @param spde A SPDE object from inla.spde2.matern() function.
#' @param Xcros A crossproduct of design matrix.
#' @param Xycros A crossproduct of design matrix and BOLD y.
#' @param thresholds A vector of thresholds for activation maps.
#' @param ind_beta A vector of indices of beta.
#' @return A list containing...
#' @export
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#' @importFrom excursions excursions.mc
#'
#' @examples \dontrun{}
beta.posterior.thetasamp <- function(theta, spde, Xcros, Xycros, thresholds, alpha=0.01, ind_beta){

	# print('Constructing joint precision')
	prec.error <- exp(theta[1])
	K <- length(theta[-1])/2
	M <- length(Xcros)

	#construct prior precision matrix for beta, Q_theta,
	#for given sampled values of theta
	theta.beta <- list()
	Q.beta <- list()
	for(k in 1:K) {
		theta.beta[[k]] <- theta[(2:3) + 2*(k-1)] #2:3, 4:5, ...
		Q.beta[[k]] <- inla.spde2.precision(spde, theta = theta.beta[[k]])
	}
	Q <- bdiag(Q.beta)

	beta.samp.pop <- 0
	beta.mean.pop <- 0
	#~25 seconds per subject
	# print('Looping over subjects')
	for(mm in 1:M){
		Xcros.mm <- Xcros[[mm]]
		Xycros.mm <- Xycros[[mm]]
		Q.m <- prec.error*Xcros.mm + Q
		mu.m <- inla.qsolve(Q.m, prec.error*Xycros.mm) #20 sec
		beta.mean.pop <- beta.mean.pop + mu.m/M

		#draw samples from pi(beta_m|theta,y)

	  	#this only works when the pop-level quantity of interest is the average activation for each task
	  	#can use the linear combination matrix A to make more general
	  	#n=10: 20 sec, n=100: 90 sec
	    beta.samp.m <- inla.qsample(n = 100, Q = Q.m, mu = mu.m) #NKx100
	    beta.samp.pop <- beta.samp.pop + beta.samp.m/M #NKx100
	}
	mu.theta <- matrix(beta.mean.pop, ncol=1)

	#3.5-7 seconds per activation threshold
	# print('Looping over activation thresholds')
	n.mesh <- spde$n.spde
	U <- length(thresholds)
	F.theta <- vector('list', U)
  	for(u in 1:U){
  		F.theta[[u]] <- matrix(nrow=n.mesh, ncol=K)
  		thr <- thresholds[u]
  		for(k in 1:K){
  			res_beta.theta.k <- excursions.mc(beta.samp.pop, u = thr, ind = ind_beta[[k]], type = '>', alpha = alpha, verbose = FALSE)
  			F.theta[[u]][,k] <- res_beta.theta.k$F[ind_beta[[k]]]
  		}
  		F.theta[[u]][is.na(F.theta[[u]])] <- 0
  	}

  	result <- list(mu.theta, F.theta)
  	names(result) <- c('mu','F')
  	return(result)
}
