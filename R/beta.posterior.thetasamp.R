#' Internal function used in joint approach to group-analysis
#'
#' @param theta A sample of theta (hyperparameters).
#' @param spde A SPDE object from inla.spde2.matern() function.
#' @param Xcros A crossproduct of design matrix.
#' @param Xycros A crossproduct of design matrix and BOLD y.
#' @param contrasts A list of vectors of length M*K specifying the contrasts of interest.  See Details for more information.
#' @param quantiles Vector of posterior quantiles to return in addition to the posterior mean
#' @param excursion_type Vector of excursion function type (">", "<", "!=") for each contrast
#' @param gamma Vector of activation thresholds for each contrast
#' @param alpha Significance level for activation for the excursion sets
#' @param nsamp_beta The number of samples to draw from p(beta|theta,y)
#'
#' @details The contrast vector specifies the group-level quantity of interest.  For example, the vector `rep(1/M,M*K)` would return the group average for each of K tasks;
#' the vector `c(rep(1/M1,M1*K)`, `rep(-1/M2,M2*K))` would return the difference between the average within two groups of size M1 and M2, respectively, for each of K tasks;
#' the vector `rep(rep(1/M,-1/M,0,...,0),each=V),M)` would return the difference between the first two tasks, averaged over all subjects.
#' @return A list containing...
#' @export
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#' @importFrom excursions excursions.mc
#' @importFrom Matrix Diagonal chol
#' @importFrom INLA inla.spde2.precision inla.qsample
#'
beta.posterior.thetasamp <- function(theta, spde, Xcros, Xycros, contrasts, quantiles, excursion_type, gamma, alpha, nsamp_beta=100){

  n.mesh <- spde$n.spde

  # print('contructing joint precision')
  prec.error <- exp(theta[1])
  K <- length(theta[-1])/2
  M <- length(Xcros)

  #contruct prior precision matrix for beta, Q_theta for given sampled value of thetas
  Q.beta <- list()
  for(k in 1:K) {
    theta_k <- theta[(2:3) + 2*(k-1)] #1:2, 2:3, 4:5, ...
    Q.beta[[k]] <- inla.spde2.precision(spde, theta = theta_k)
  }
  Q <- bdiag(Q.beta)

  # #make block diagonal Q in multi-session case
  # tmp1 <- nrow(Q)
  # tmp2 <- nrow(Xcros[[1]])
  # n_sess <- tmp2/tmp1
  # if(n_sess != round(n_sess)) stop('Something is wrong. Dimensions of prior precision and X matrix not compatible.')
  # Q <- bdiag(rep(list(Q), n_sess))
  #
  N <- dim(Q.beta[[1]])[1]
  if(N != n.mesh) stop('Length of betas does not match number of vertices in mesh. Contact developer to request this functionality.')

  #beta.samp.pop <- 0
  #beta.mean.pop <- 0
  #beta.mean.pop.mat <- NULL
  beta.samples <- NULL
  #~5 seconds per subject with PARDISO
  # print('Looping over subjects or sessions')
  for(mm in 1:M){

    #compute posterior mean and precision of beta|theta
    Xcros.mm <- Xcros[[mm]] #X'X
    Xycros.mm <- Xycros[[mm]] #X'y
    Q.m <- prec.error*Xcros.mm + Q
    mu.m <- inla.qsolve(Q.m, prec.error*Xycros.mm) #NK x 1 -- 2 minutes, but only 2 sec with PARDISO!

    #draw samples from pi(beta_m|theta,y)
    beta.samp.m <- inla.qsample(n = nsamp_beta, Q = Q.m, mu = mu.m) #NK x nsamp_beta  -- 2 minutes, but only 2 sec with PARDISO!

    #concatenate samples over models
    beta.samples <- rbind(beta.samples, beta.samp.m) #(N*K*M) x nsamp_beta

    # ## group mean
    # beta.mean.pop <- beta.mean.pop + mu.m/M
    # beta.samp.pop <- beta.samp.pop + beta.samp.m/M #NKx100

    #       ## contrast
    # 	    beta.mean.pop.mat <- c(beta.mean.pop.mat, mu.m) #NKM x 1
    # 	    beta.samp.pop.mat <- rbind(beta.samp.pop.mat, beta.samp.m) #NKM x nsamp_beta

  }

  ## Compute results for beta averages over subjects (default)
  #mu.theta <- matrix(beta.mean.pop, ncol=1)

  # 	if(excursion_type != 'none'){
  #
  #   	U <- length(gamma)
  #   	F.theta <- vector('list', U)
  #   	for(u in 1:U){
  #   		F.theta[[u]] <- matrix(nrow=n.mesh, ncol=K)
  #   		thr <- gamma[u]
  #   		for(k in 1:K){
  #   			res_beta.theta.k <- excursions.mc(beta.samp.pop, u = thr, ind = ind_beta[[k]], type = excursion_type, alpha = alpha)
  #   			F.theta[[u]][,k] <- res_beta.theta.k$F[ind_beta[[k]]]
  #   		}
  #   		F.theta[[u]][is.na(F.theta[[u]])] <- 0
  #   	}
  # 	}

  if(excursion_type[1] == 'none') do_excur <- FALSE else do_excur <- TRUE

  # Loop over contrasts
  num_contrasts <- length(contrasts)
  mu.contr <- F.contr <- matrix(NA, nrow=n.mesh, ncol=num_contrasts)
  num_quantiles <- length(quantiles)
  quantiles.contr <- rep(list(matrix(NA, nrow=n.mesh, ncol=num_contrasts)), num_quantiles)
  names(quantiles.contr) <- paste0('q',quantiles)
  for(l in 1:num_contrasts){

    #Construct "A" matrix from paper
    ctr.mat <- kronecker(t(contrasts[[l]]), Diagonal(n.mesh, 1))

    #beta.mean.pop.contr <- as.vector(ctr.mat%*%beta.mean.pop.mat)  # NKx1 or Nx1
    samples_l <- as.matrix(ctr.mat%*%beta.samples)  # NKx100 or Nx100
    mu.contr[,l] <- rowMeans(samples_l)
    if(num_quantiles > 0){
      for(iq in 1:num_quantiles){
        quantiles.contr[[iq]][,l] <- apply(samples_l, 1, quantile, quantiles[iq])
      }
    }

    # Estimate excursions set for current contrast
    if(do_excur){
      type_l <- excursion_type[l]
      thr_l <- gamma[l]
      excur_l <- excursions.mc(samples_l, u = thr_l, type = type_l, alpha = alpha)
      F.contr[,l] <- excur_l$F
    }
  }

  # 	result <- list(mean = list(mu = mu.theta, F = F.theta),
  #                 contr = list(mu = mu.contr, F = F.contr))
  result <- list(mu = mu.contr, quantiles=quantiles.contr, F = F.contr)
  return(result)
}
