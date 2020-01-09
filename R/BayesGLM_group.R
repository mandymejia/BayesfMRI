#' Applies joint approach to group-level analysis to task fMRI data
#'
#' @param result Either (1) a list of length M of objects of class BayesGLM, or (2) a character vector of length M of file names output from the BayesGLM function. M is the number of subjects.
#' @param A A matrix to translate between original and mesh locations.
#' @param contrasts A vector of length M*K specifying the contrast of interest.  See Details for more information.
#' @param no_cores The number of cores to use for sampling in parallel
#' @param type The type of excursion function for mean beta (">", "<", "!=")
#' @param type.contr The list of types of excursion function for contrasts (">", "<", "!=")
#' @param thresholds The vector of activation thresholds
#' @param thresholds.contr The list of activation thresholds for contrasts
#' @param alpha The significance level for activation
#' @param alpha.contr The list of activation significance levels for contrasts
#'
#' @details The contrast vector specifies the group-level quantity of interest.  For example, the vector `rep(1/M,M*K)` would return the group average for each of K tasks;
#' the vector `c(rep(1/M1,M1*K)`, `rep(-1/M2,M2*K))` would return the difference between the average within two groups of size M1 and M2, respectively, for each of K tasks;
#' the vector `rep(rep(1/M,-1/M,0,...,0)),M)` would return the difference between the first two tasks, averaged over all subjects.
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern
#' @importFrom MASS mvrnorm
#' @import parallel
#'
#' @examples \dontrun{}
BayesGLM_group <- function(result, A, contrasts = NULL, thresholds = 0, thresholds.contr = NULL, type = NULL, type.contr = NULL, alpha = 0.05, alpha.contr = NULL, no_cores=NULL){

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
    beta.post.samps <- apply(theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, Xcros.all, Xycros.all, contrasts=contrasts, thresholds=thresholds, thresholds.contr=thresholds.contr, type=type, type.contr=type.contr, alpha=alpha, alpha.contr=alpha.contr, ind_beta=ind_beta)
  } else {
    ## Not sure how to change the following
    max_no_cores <- min(detectCores() - 1, 25)
    no_cores <- min(max_no_cores, no_cores)
    cl <- makeCluster(no_cores)
    beta.post.samps <- parApply(cl, theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, K=K, M, Xcros.all, Xycros.all, thresholds=thresholds, alpha=alpha, ind_beta=ind_beta)
    print(Sys.time() - t0)
    stopCluster(cl)
  }


### Computing posterior quantities of beta averaged over subjects (summing over theta)
#organize samples
U <- length(thresholds)
mu.tot <- matrix(nrow=K*n.mesh, ncol=nsamp)
F.tot <- rep(list(rep(list(matrix(nrow=n.mesh, ncol=nsamp)), K)), U) #for each activation threshold and task, a Vx50 matrix
for(itheta in 1:nsamp){
  mu.tot[,itheta] <- beta.post.samps[[itheta]]$mean$mu
  for(u in 1:U){
    for(k in 1:K){
      F.tot[[u]][[k]][,itheta] <- beta.post.samps[[itheta]]$mean$F[[u]][,k]
    }
  }
}
## Sum over samples using weights
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



### Computing posterior quantities of contrasts (summing over theta)
if(is.null(contrasts) == FALSE){
#organize samples
betas.all.contr <- probs.all.contr <- active.all.contr <- list()
 for(i.contr in 1:length(contrasts)){
  K.contr <- dim(contrasts[[i.contr]])[1]/n.mesh
  U.contr <- length(thresholds.contr[[i.contr]])
  mu.tot.contr <- matrix(nrow=K.contr*n.mesh, ncol=nsamp)
  F.tot.contr <- rep(list(rep(list(matrix(nrow=n.mesh, ncol=nsamp)), K.contr)), U.contr) #for each activation threshold and contrast, a Vx50 matrix
  
  for(itheta in 1:nsamp){
    mu.tot.contr[,itheta] <- beta.post.samps[[itheta]]$contr$mu[[i.contr]]
    for(u in 1:U.contr){
      for(k in 1:K.contr){
        F.tot.contr[[u]][[k]][,itheta] <- beta.post.samps[[itheta]]$contr$F[[i.contr]][[u]][,k]
      }
    }
  }
  ## Sum over samples using weights
  betas.all.c <- matrix(0, nrow=n.mesh, ncol=K.contr)
  probs.all.c <- active.all <- array(0, dim=c(n.mesh, K.contr, U.contr)) #last dimension is for different activation thresholds
  
  #posterior mean
  beta.pop.contr <- as.vector(mu.tot.contr%*%wt)
  for(k in 1:K.contr){
    beta.pop.k <- beta.pop.contr[ind_beta[[k]]]
    betas.all.c[,k] <- as.vector(beta.pop.k)
  }
  #posterior probabilities and excursion sets
  for(u in 1:U.contr){
    for(k in 1:K.contr){
      F.pop.uk <- as.vector(F.tot.contr[[u]][[k]]%*%wt)
      E.pop.uk <- rep(0, length(F.pop.uk))
      E.pop.uk[F.pop.uk > 1 - alpha.contr[i.contr]] <- 1
      probs.all.c[,k,u] <- as.vector(F.pop.uk)
      active.all.c[,k,u] <- as.vector(E.pop.uk)
    }
  }
  betas.all.contr[[i.contr]] <- betas.all.c
  probs.all.contr[[i.contr]] <- probs.all.c
  active.all.contr[[i.contr]] <- active.all.c
 }
} else{
  betas.all.contr <- probs.all.contr <- active.all.contr <- NULL
}

### Save all results 
result <- list(mean = list(beta_estimates = betas.all, ppm = probs.all, active = active.all), contrasts = list(beta_estimates = betas.all.contr, ppm = probs.all.contr, active = active.all.contr))

class(result) <- "BayesGLM_group"

return(result)

}


