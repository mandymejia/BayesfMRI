#' Internal function used in joint approach to group-analysis
#'
#' @param theta A sample of theta (hyperparameters)
#' @param spde A SPDE object from inla.spde2.matern() function.
#' @param Xcros A crossproduct of design matrix.
#' @param Xycros A crossproduct of design matrix and BOLD y.
#' @param thresholds A vector of thresholds for activation maps.
#' @param ind_beta A vector of indices of beta.
#' @param contrasts A list of vectors of length M*K specifying the contrasts of interest.  See Details for more information.
#' @details The contrast vector specifies the group-level quantity of interest.  For example, the vector `rep(1/M,M*K)` would return the group average for each of K tasks;
#' the vector `c(rep(1/M1,M1*K)`, `rep(-1/M2,M2*K))` would return the difference between the average within two groups of size M1 and M2, respectively, for each of K tasks;
#' the vector `rep(rep(1/M,-1/M,0,...,0),each=V),M)` would return the difference between the first two tasks, averaged over all subjects.
#' @return A list containing...
#' @export
#' @importFrom excursions excursions.mc
#'
#' @examples \dontrun{}
beta.posterior.thetasamp <- function(theta, spde, Xcros, Xycros, contrasts, thresholds, thresholds.contr, type, type.contr, alpha, alpha.contr, ind_beta){

	# print('contructing joint precision')
	prec.error <- exp(theta[1])
	K <- length(theta[-1])/2
	M <- length(Xcros)

	#contruct prior precision matrix for beta, Q_theta,
	#for given sampled values of theta
	theta.beta <- list()
	Q.beta <- list()
	for(k in 1:K) {
		theta.beta[[k]] <- theta[(2:3) + 2*(k-1)] #2:3, 4:5, ...
		Q.beta[[k]] <- inla.spde2.precision(spde, theta = theta.beta[[k]])
	}
	Q <- bdiag(Q.beta)
	N <- dim(Q.beta[[1]])[1]
	Idn <- Diagonal(N, x = 1)

	beta.samp.pop <- 0
	beta.mean.pop <- 0
	beta.mean.pop.mat <- NULL
	beta.samp.pop.mat <- NULL
	#~25 seconds per subject
	# print('Looping over subjects')
	for(mm in 1:M){
		Xcros.mm <- Xcros[[mm]]
		Xycros.mm <- Xycros[[mm]]
		Q.m <- prec.error*Xcros.mm + Q
		mu.m <- inla.qsolve(Q.m, prec.error*Xycros.mm) #20 sec  NK x 1
		#draw samples from pi(beta_m|theta,y)
	  beta.samp.m <- inla.qsample(n = 100, Q = Q.m, mu = mu.m) #NK x 100
	  
	  if(is.null(contrasts)){ ## return group mean for each task by default
	    beta.mean.pop <- beta.mean.pop + mu.m/M
	    beta.samp.pop <- beta.samp.pop + beta.samp.m/M #NKx100
	  } else{
	    beta.mean.pop.mat <- c(beta.mean.pop.mat, mu.m) #NKM x 1
	    beta.samp.pop.mat <- rBind(beta.samp.pop.mat, beta.samp.m) #NKM x 100
	  }
	}
	
	## Compute results for beta averages over subjects (default)
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
  			res_beta.theta.k <- excursions.mc(beta.samp.pop, u = thr, ind = ind_beta[[k]], type = type, alpha = alpha, verbose = FALSE)
  			F.theta[[u]][,k] <- res_beta.theta.k$F[ind_beta[[k]]]
  		}
  		F.theta[[u]][is.na(F.theta[[u]])] <- 0
  	}
	
	
	## Compute results for contrasts
	if(is.null(contrasts) == FALSE){
	  mu.contr <- F.contr <- vector('list', length(contrasts))
	  
	  for(n.ctr in 1:length(contrasts)){
	    
	    # Loop over contrasts
	    ctr.vec <- contrasts[[n.ctr]]
	    beta.mean.pop.contr <- as.vector(t(kronecker(ctr.vec, Idn))%*%beta.mean.pop.mat)  # NKx1 or Nx1
	    beta.samp.pop.contr <- t(kronecker(ctr.vec, Idn))%*%beta.samp.pop.mat  # NKx100 or Nx100
	    mu.theta.contr <- matrix(beta.mean.pop.contr, ncol=1)
	    
	    # Loop over activation thresholds for each contrast
	    n.mesh <- spde$n.spde
	    K.contr <- dim(mu.theta.contr)[1]/n.mesh
	    U <- length(thresholds.contr)
	    F.theta.contr <- vector('list', U)
	    for(u in 1:U){
	      F.theta.contr[[u]] <- matrix(nrow = n.mesh, ncol = K.contr)
	      thr <- thresholds.contr[u]
	      for(k in 1:K.contr){
	        res_contr.k <- excursions.mc(beta.samp.pop.contr, u = thr, ind = ind_beta[[k]], type = type.contr, alpha = alpha.contr, verbose = FALSE)
	        F.theta.contr[[u]][,k] <- res_contr.k$F[ind_beta[[k]]]
	      }
	      F.theta.contr[[u]][is.na(F.theta.contr[[u]])] <- 0
	    }
	    
	    mu.contr[[n.ctr]] <- mu.theta.contr
	    F.contr[[n.ctr]] <- F.theta.contr
	  }
	  
	}else{
	  mu.contr = NULL
	  F.contr = NULL
	}

  	result <- list(mean = list(mu = mu.theta, F = F.theta), contr = list(mu = mu.contr, F = F.contr))
  	# names(result) <- c('mu','F')
  	return(result)
}
