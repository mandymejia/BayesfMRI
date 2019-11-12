#' Applies joint approach to group-level analysis to task fMRI data
#'
#' @param results Either (1) a list of length M of objects of class BayesGLM, or (2) a character vector of length M of file names output from the BayesGLM function. M is the number of subjects.
#' @param A A matrix to translate between original and mesh locations.
#' @param contrast A vector specifying the contrast of interest.  See Details for more information.
#' @param no_cores The number of cores to use for sampling in parallel
#' @param thresholds The vector of activation thresholds
#' @param alpha The significance level for activation
#'
#' @details The contrast vector specifies the group-level quantity of interest.  For example, the vector rep(1,M*K) would return the group average for each of K tasks;
#' the vector `c(rep(1,M1*K)`, `rep(-1,M2*K))` would return the difference between the average within two groups of size M1 and M2, respectively, for each of K tasks;
#' the vector `rep(rep(1,-1,0,0),each=V),M)` would return the difference between the first two tasks (of 4), averaged over all subjects.
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern
#' @importFrom MASS mvrnorm
#' @import parallel
#'
#' @examples \dontrun{}
BayesGLM_group <- function(result, A, contrasts = NULL, thresholds = 0, alpha = 0.05, no_cores=NULL){

  require(INLA)
  # Find the numnber of subjects.
  subject_names <- names(result)
  M <- length(subject_names)

  spde <- inla.spde2.matern(result[[1]]$mesh)

  # Collecting theta posteriors from subject models
  theta.sub <- NULL
  mu.theta.tmp <- Q.theta <- 0

  for(m in 1:M){

      # Save it in BayesGLM()
      res.hyper <- result[[m]]$INLA_result$summary.hyperpar

      mu.tmp <- result[[m]]$mu.theta
      Q.tmp <- result[[m]]$Q.theta

      mu.theta.tmp <- mu.theta.tmp + as.vector(Q.tmp%*%mu.tmp)
      Q.theta <- Q.theta + Q.tmp
      theta.sub <- cbind(theta.sub, res.hyper$mode)
      rm(mu.tmp, Q.tmp)
    }
    mu.theta <- solve(Q.theta, mu.theta.tmp)

  # Drawing samples from q(theta|y)
  nsamp <- 50
  logwt <- rep(NA, nsamp)
  theta.tmp <- MASS::mvrnorm(nsamp, mu.theta, solve(Q.theta))
  for(i in 1:nsamp){
    logwt[i] <- F.logwt(theta.tmp[i,], spde, mu.theta, Q.theta, M)
  }

  #weights to apply to each posterior sample of theta
  wt.tmp <- exp(logwt - max(logwt))
  wt <- wt.tmp/(sum(wt.tmp))


  ## Create index vectors
  K <- (dim(theta.tmp)[2] - 1)/2
  n.mesh <- mesh$n
  ind_beta <- list()
  for(k in 1:K){
    ind_beta[[k]] <- 1:n.mesh + (k-1)*n.mesh
  }

  ## Compute cross-products for single session
  A.lst <- vector("list", K)
  for(k in 1:K){
    A.lst[[k]] <- A
  }
  Amat.tot <- bdiag(A.lst)

  Xcros.all <- Xycros.all <- vector("list", M)
  for(m in 1:M){
    y_vec <- result[[m]]$y
    X_list <- result[[m]]$X

    Xmat <- X_list[[1]]%*%Amat.tot
    Xcros.all[[m]] <- Matrix::crossprod(Xmat)
    Xycros.all[[m]] <- Matrix::crossprod(Xmat, y_vec)
  }

  #get posterior quantities of beta, conditional on a value of theta
  if(is.null(no_cores)) {
    parallel <- FALSE
    beta.post.samps <- apply(theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, Xcros.all, Xycros.all, thresholds=thresholds, alpha=alpha, ind_beta=ind_beta)
  } else {
    max_no_cores <- min(detectCores() - 1, 25)
    no_cores <- min(max_no_cores, no_cores)
    cl <- makeCluster(no_cores)
    beta.post.samps <- parApply(cl, theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, K=K, M, Xcros.all, Xycros.all, thresholds=thresholds, alpha=alpha, ind_beta=ind_beta)
    print(Sys.time() - t0)
    stopCluster(cl)
  }

  #organize samples
  U <- length(thresholds)
  mu.tot <- matrix(nrow=K*n.mesh, ncol=nsamp)
  F.tot <- rep(list(rep(list(matrix(nrow=n.mesh, ncol=nsamp)), K)), U) #for each activation threshold and task, a Vx50 matrix
  for(itheta in 1:nsamp){
    mu.tot[,itheta] <- beta.post.samps[[itheta]]$mu
    for(u in 1:U){
      for(k in 1:K){
        F.tot[[u]][[k]][,itheta] <- beta.post.samps[[itheta]]$F[[u]][,k]
      }
    }
  }

  # Computing posterior quantities of beta, summing over theta')

  ### Sum over samples using weights, combine hemispheres (< 1 sec)
  betas.all <- matrix(0, nrow=n.mesh, ncol=K)
  probs.all <- active.all <- array(0, dim=c(n.mesh, K, U)) #last dimension is for different activation thresholds

  #posterior mean
  beta.pop <- as.vector(mu.tot%*%wt)

  for(k in 1:K){
      beta.pop.k <- beta.pop[ind_beta[[k]]]
      betas.all[,k] <- as.vector(beta.pop.k)
  }

  #posterior probabilities and excursion sets
  for(u in 1:U){
    for(k in 1:K){
      F.pop.uk <- as.vector(F.tot[[u]][[k]]%*%wt)
      E.pop.uk <- rep(0, length(F.pop.uk))
      E.pop.uk[F.pop.uk > 1-alpha] <- 1
      probs.all[,k,u] <- as.vector(F.pop.uk)
      active.all[,k,u] <- as.vector(E.pop.uk)
    }
  }


  result <- list(beta_estimates = betas.all, ppm = probs.all, active = active.all)

  class(result) <- "BayesGLM_group"

  return(result)

}


