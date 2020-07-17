#' Computes Bayesian GLM posterior means, quantiles and excursions sets (areas of activation) for summaries or contrasts across sessions within a single multi-session model
#'
#' @param result An object of class BayesGLM containing the results of a multi-session analysis. Must contain INLA_result!  (See `help(BayesGLM)`.)
#' @param contrasts A list of vectors, each length M*K, specifying the contrast(s) of interest across sessions, where M is the number of sessions analyzed and K is the number of tasks.  See Details for more information.
#' @param quantiles Vector of posterior quantiles to return in addition to the posterior mean
#' @param excursion_type The type of excursion function for the contrast (">", "<", "!="), or a vector thereof (each element corresponding to one contrast).  If none provided, no inference performed.
#' @param gamma Activation threshold for the excursion set, or a vector thereof (each element corresponding to one contrast).
#' @param alpha Significance level for activation for the excursion set
#' @param no_cores The number of cores to use for sampling in parallel. If NULL, do not run in parallel.
#' @param nsamp_theta Number of theta values to sample from posterior
#' @param nsamp_beta Number of beta vectors to sample conditional on each theta value sampled
#'
#' @details Performs posterior inference on summaries or contrasts across the sessions in a multi-session Bayesian GLM.  Let M be the number of sessions analyzed in the
#' multi-session Bayesian GLM, and let K be the number of tasks.  Computes the posterior mean of each contrast and produces excursions sets representing areas of
#' activation (above a specified activation threshold gamma) for each contrast map based on the joint (across locations) posterior distribution.
#'
#' Each contrast vector specifies an across-session summary of interest.  For example, the contrast vector `rep(c(1/M, rep(0, M-1)), K)` represents the across-session average
#' for the first task; `c( rep(c(1/M1, rep(0, M1-1)), K), rep(c(-1/M2, rep(0, M2-1)), K))` represents the difference between the first M1 sessions and the remaining M2 sessions (M1+M2=M) for the first task;
#'`rep( c(1/M,-1/M,rep(0, K-2)) ,M)` represents the difference between the first two tasks, averaged over all sessions.
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern inla.pardiso.check inla.setOption inla.qsample
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag crossprod
#' @importFrom abind abind
#' @import parallel
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
BayesGLM_group.within <- function(result,
                                     contrasts = NULL,
                                     quantiles = NULL,
                                     excursion_type,
                                     gamma = 0,
                                     alpha = 0.05,
                                     no_cores = NULL,
                                     nsamp_theta = 50,
                                     nsamp_beta = 100){

  flag <- inla.pardiso.check()
  if(grepl('FAILURE',flag)) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO is required for computational efficiency. See inla.pardiso().')
  inla.setOption(smtp='pardiso')

  # Find the numnber of sessions and tasks
  session_names <- names(result$design)
  M <- length(session_names)
  K <- ncol(result$design[[1]])

  #### SET UP OR CHECK CONTRAST VECTOR(S)

  # If contrasts is null, by default set up a contrast vector that will compute the average across subjects/sessions for each task
  if(!is.null(contrasts) & !is.list(contrasts)) contrasts <- list(contrasts)
  if(is.null(contrasts)) {
    contrasts <- vector('list', length=K)
    for(k in 1:K){
      contrast_onesubj_k <- c(rep(0, k-1), 1/M, rep(0, K-k)) #(1/M, 0, ... 0) for k=1, (0, 1/M, 0, ..., 0) for k=2, ..., (0, ..., 0, 1/M) for k=K
      contrast_allsubj_k <- rep(contrast_onesubj_k, M)
      contrasts[[k]] <- contrast_allsubj_k
    }
  }

  #Check that each contrast vector is numeric and length M*K
  num_contrasts <- length(contrasts)
  length_each_contrast <- sapply(contrasts, length)
  class_each_contrast <- sapply(contrasts, is.numeric)
  if(any(length_each_contrast != M*K)) stop('each contrast vector must be of length M*K')
  if(any(!class_each_contrast)) stop('each contrast vector must be numeric, but at least one is not')

  #Check excursions settings
  if(exists(excursion_type)) {
    do_excur <- TRUE
    if(length(excursion_type) == 1) excursion_type <- rep(excursion_type, num_contrasts)
    if(length(gamma) == 1) gamma <- rep(gamma, num_contrasts)
    if(length(alpha) != 1) stop('alpha must be a scalar value')
    if(length(gamma) != num_contrasts) stop('Length of gamma must match number of contrasts or be equal to one.')
    if(length(excursion_type) != num_contrasts) stop('Length of excursion_type must match number of contrasts or be equal to one.')
  } else {
    excursion_type <- 'none'
  }

  #Check quantiles argument
  if(!is.null(quantiles)){
    if(any(quantiles > 1 | quantiles < 0)) stop('All elements of quantiles argument must be between 0 and 1.')
  }

  #### LOOP OVER HEMISPHERES

  do_left <- do_right <- TRUE
  if(is.null(result$GLMs_Bayesian$cortex_left)) do_left <- FALSE
  if(is.null(result$GLMs_Bayesian$cortex_right)) do_right <- FALSE

  #put each hemisphere results into a list to create cifti objects later
  posterior_means <- vector('list', length=2)
  posterior_quantiles <- vector('list', length=2)
  posterior_probs <- vector('list', length=2)
  posterior_active <- vector('list', length=2)

  for(h in 1:2){

    if(h==1){
      if(!do_left) next
      result_h <- result$GLMs_Bayesian$cortex_left
    }
    if(h==1){
      if(!do_right) next
      result_h <- result$GLMs_Bayesian$cortex_right
    }

    # Mesh and SPDE object
    mesh <- result_h$mesh
    spde <- inla.spde2.matern(mesh)
    mask <- result_h$mask

    #### GET THETA POSTERIOR FROM MULTI-SESSION MODEL

    Q.theta <- result_h$Q.theta
    mu.theta <- result_h$mu.theta

    #### DRAW SAMPLES FROM POSTERIOR OF THETA (Normal approximation)

    theta.samp <- inla.qsample(n=nsamp_theta, Q = Q.theta, mu = mu.theta)

    # ## Create index vectors
    # K <- (dim(theta.samp)[2] - 1)/2
    # n.mesh <- mesh$n
    # ind_beta <- list()
    # for(k in 1:K){
    #   ind_beta[[k]] <- 1:n.mesh + (k-1)*n.mesh
    # }

    # ## Compute cross-products for single session
    # A.lst <- vector("list", K)
    # for(k in 1:K){
    #   A.lst[[k]] <- A
    # }
    # Amat.tot <- bdiag(A.lst)

    Xcros <- Xycros <- vector('list', M)
    y_vec <- result_h$y
    X_list <- result_h$X #%*%Amat.tot #for multi-session data, make X a block diagonal matrix
    for(mm in 1:M){
      inds_m <- (1:nrow(X_list[[1]])) + (mm-1)*nrow(X_list[[1]])
      Xcros[[mm]] <- crossprod(X_list[[mm]])
      Xycros[[mm]] <- crossprod(X_list[[mm]], y_vec[inds_m])
    }

    print('Computing posterior quantities of beta for each value of theta')

    #get posterior quantities of beta, conditional on a value of theta
    #if(is.null(no_cores)) {

      #apply is crashing R!  Try looping?
      beta.posteriors <- vector('list', length=nsamp_theta)
      for(isamp in 1:nsamp_theta){
        print(paste0('theta value ',isamp,' of ',nsamp_theta))
        print(
          system.time(
            beta.posteriors[[isamp]] <- beta.posterior.thetasamp(theta.samp[,isamp],
                                 spde=spde,
                                 Xcros=Xcros,
                                 Xycros=Xycros,
                                 contrasts=contrasts,
                                 quantiles=quantiles,
                                 excursion_type=excursion_type,
                                 gamma=gamma,
                                 alpha=alpha,
                                 nsamp_beta=nsamp_beta))) #20 sec/iter
      }

    #   beta.posteriors <- apply(theta.samp, MARGIN=2,
    #                            FUN=beta.posterior.thetasamp,
    #                            spde=spde,
    #                            Xcros=Xcros,
    #                            Xycros=Xycros,
    #                            contrasts=contrasts,
    #                            quantiles=quantiles,
    #                            excursion_type=excursion_type,
    #                            gamma=gamma,
    #                            alpha=alpha)
    #
    # } else {
    #   max_no_cores <- min(detectCores() - 1, 25)
    #   no_cores <- min(max_no_cores, no_cores)
    #   cl <- makeCluster(no_cores)
    #   beta.posteriors <- parApply(cl, theta.samp, MARGIN=2,
    #                            FUN=beta.posterior.thetasamp,
    #                            spde=spde,
    #                            Xcros=Xcros,
    #                            Xycros=Xycros,
    #                            contrasts=contrasts,
    #                            quantiles=quantiles,
    #                            excursion_type=excursion_type,
    #                            gamma=gamma,
    #                            alpha=alpha)
    #   stopCluster(cl)
    # }


    # ### Computing posterior quantities of beta averaged over subjects (summing over theta)
    # #organize samples
    # mu.tot <- matrix(nrow=K*n.mesh, ncol=nsamp_theta)
    # F.tot <- rep(list(rep(list(matrix(nrow=n.mesh, ncol=nsamp_theta)), K)), U) #for each activation threshold and task, a Vx50 matrix
    # for(itheta in 1:nsamp_theta){
    #   mu.tot[,itheta] <- beta.posteriors[[itheta]]$mean$mu
    #   for(u in 1:U){
    #     for(k in 1:K){
    #       F.tot[[u]][[k]][,itheta] <- beta.posteriors[[itheta]]$mean$F[[u]][,k]
    #     }
    #   }
    # }
    # ## Sum over samples using weights
    # betas.all <- matrix(0, nrow=n.mesh, ncol=K)
    # probs.all <- active.all <- array(0, dim=c(n.mesh, K, U)) #last dimension is for different activation thresholds

    #posterior mean (average over all thetas)
    mu.bytheta <- lapply(beta.posteriors, function(x) return(x$mu)) #extract conditional posterior means
    mu.bytheta <- abind(mu.bytheta, along=3) #list to array
    mu.all <- matrix(NA, nrow=length(mask), ncol = num_contrasts)
    mu.all[mask,] <- apply(mu.bytheta, c(1,2), mean) #average over thetas


    # #posterior mean
    # beta.pop <- as.vector(mu.tot%*%wt)
    # for(k in 1:K){
    #   beta.pop.k <- beta.pop[ind_beta[[k]]]
    #   betas.all[,k] <- as.vector(beta.pop.k)
    # }

    #posterior quantiles
    num_quantiles <- length(quantiles)
    if(num_quantiles > 0){
      quantiles.all <- vector('list', length=num_quantiles)
      names(quantiles.all) <- paste0('q',quantiles)
      for(iq in 1:num_quantiles){
        quantiles.bytheta <- lapply(beta.posteriors, function(x) return(x$quantiles[[iq]]))
        quantiles.bytheta <- abind(quantiles.bytheta, along=3)
        quantiles.all[[iq]] <- matrix(NA, nrow=length(mask), ncol = num_contrasts)
        quantiles.all[[iq]][mask,] <- apply(quantiles.bytheta, c(1,2), mean)
      }
    }

    #posterior probabilities and excursion sets
    F.bytheta <- lapply(beta.posteriors, function(x) return(x$F))
    F.bytheta <- abind(F.bytheta, along=3)
    F.all <- matrix(NA, nrow=length(mask), ncol = num_contrasts)
    F.all[mask,] <- apply(F.bytheta, c(1,2), mean) #joint PPM / posterior probabilities of activation
    E.all <- F.all*0
    E.all[F.all > 1-alpha] <- 1 #active areas
    active.all <- matrix(NA, nrow=length(mask), ncol = num_contrasts)
    active.all[mask,] <- E.all


    # for(u in 1:U){
    #   for(k in 1:K){
    #     F.pop.uk <- as.vector(F.tot[[u]][[k]]%*%wt)
    #     E.pop.uk <- rep(0, length(F.pop.uk))
    #     E.pop.uk[F.pop.uk > 1-alpha] <- 1
    #     probs.all[,k,u] <- as.vector(F.pop.uk)
    #     active.all[,k,u] <- as.vector(E.pop.uk)
    #   }
    # }
    #


    # ### Computing posterior quantities of contrasts (summing over theta)
    # if(is.null(contrasts) == FALSE){
    #   #organize samples
    #   betas.all.contr <- probs.all.contr <- active.all.contr <- list()
    #   for(i.contr in 1:length(contrasts)){
    #     K.contr <- dim(contrasts[[i.contr]])[1]/n.mesh
    #     U.contr <- length(thresholds.contr[[i.contr]])
    #     mu.tot.contr <- matrix(nrow=K.contr*n.mesh, ncol=nsamp_theta)
    #     F.tot.contr <- rep(list(rep(list(matrix(nrow=n.mesh, ncol=nsamp_theta)), K.contr)), U.contr) #for each activation threshold and contrast, a (V x nsamp_theta) matrix
    #
    #     for(itheta in 1:nsamp_theta){
    #       mu.tot.contr[,itheta] <- beta.posteriors[[itheta]]$contr$mu[[i.contr]]
    #       for(u in 1:U.contr){
    #         for(k in 1:K.contr){
    #           F.tot.contr[[u]][[k]][,itheta] <- beta.posteriors[[itheta]]$contr$F[[i.contr]][[u]][,k]
    #         }
    #       }
    #     }
    #     ## Sum over samples using weights
    #     betas.all.c <- matrix(0, nrow=n.mesh, ncol=K.contr)
    #     probs.all.c <- active.all.c <- array(0, dim=c(n.mesh, K.contr, U.contr)) #last dimension is for different activation thresholds
    #
    #     #posterior mean
    #     beta.pop.contr <- as.vector(mu.tot.contr%*%wt)
    #     for(k in 1:K.contr){
    #       beta.pop.k <- beta.pop.contr[ind_beta[[k]]]
    #       betas.all.c[,k] <- as.vector(beta.pop.k)
    #     }
    #     #posterior probabilities and excursion sets
    #     for(u in 1:U.contr){
    #       for(k in 1:K.contr){
    #         F.pop.uk <- as.vector(F.tot.contr[[u]][[k]]%*%wt)
    #         E.pop.uk <- rep(0, length(F.pop.uk))
    #         E.pop.uk[F.pop.uk > 1 - alpha.contr[i.contr]] <- 1
    #         probs.all.c[,k,u] <- as.vector(F.pop.uk)
    #         active.all.c[,k,u] <- as.vector(E.pop.uk)
    #       }
    #     }
    #     betas.all.contr[[i.contr]] <- betas.all.c
    #     probs.all.contr[[i.contr]] <- probs.all.c
    #     active.all.contr[[i.contr]] <- active.all.c
    #   }
    # } else{
    #   betas.all.contr <- probs.all.contr <- active.all.contr <- NULL
    # }


    posterior_means[[h]] <- mu.all
    posterior_quantiles[[h]] <- quantiles.all
    posterior_probs[[h]] <- F.all
    posterior_active[[h]] <- active.all

  } #end loop over hemispheres

  ## Make cifti objects

  posterior_means_cifti <- cifti_make(cortex_left = posterior_means[[1]],
                                      cortex_right = posterior_means[[2]])
                                      #surf_left = result_svh$betas_Bayesian$LR$SURF_LEFT$surface)
  if(do_excur){
    posterior_probs_cifti <- cifti_make(cortex_left = posterior_probs[[1]],
                                        cortex_right = posterior_probs[[2]])
    posterior_active_cifti <- cifti_make(cortex_left = posterior_active[[1]],
                                         cortex_right = posterior_active[[2]])
                                         #surf_left = result_svh$betas_Bayesian$LR$SURF_LEFT$surface)
  }
  if(num_quantiles > 0){
    posterior_quantiles_cifti <- vector('list', num_quantiles)
    names(posterior_quantiles_cifti) <- paste0('q',quantiles)
    for(iq in 1:num_quantiles){
      posterior_quantiles_cifti[[1]] <- cifti_make(cortex_left = posterior_quantiles[[1]][[iq]],
                                                   cortex_right = posterior_quantiles[[2]][[iq]])
    }
  }


  ### Save all results
  result <- list(mean = posterior_means_cifti,
                 quantiles = posterior_quantiles_cifti,
                 PPM = posterior_probs_cifti,
                 active = posterior_active_cifti)

  class(result) <- "BayesGLM_group"

  return(result)

}


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
  #' @importFrom Matrix bdiag crossprod
  #' @import parallel
  #'
  #' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
  #' @examples \dontrun{}
  # BayesGLM_group <- function(result, A, contrasts = NULL, thresholds = 0, thresholds.contr = NULL, type = NULL, type.contr = NULL, alpha = 0.05, alpha.contr = NULL, no_cores=NULL){
  #
  #   # Find the numnber of subjects.
  #   subject_names <- names(result)
  #   M <- length(subject_names)
  #
  #   # Mesh and SPDE object
  #   mesh <- result[[1]]$mesh
  #   spde <- inla.spde2.matern(mesh)
  #
  #   # Collecting theta posteriors from subject models
  #   mu.theta.tmp <- Q.theta <- 0
  #
  #   for(m in 1:M){
  #
  #       mu.tmp <- result[[m]]$mu.theta
  #       Q.tmp <- result[[m]]$Q.theta
  #
  #       mu.theta.tmp <- mu.theta.tmp + as.vector(Q.tmp%*%mu.tmp)
  #       Q.theta <- Q.theta + Q.tmp
  #       rm(mu.tmp, Q.tmp)
  #     }
  #     mu.theta <- solve(Q.theta, mu.theta.tmp)
  #
  #   # Drawing samples from q(theta|y)
  #   nsamp <- 50
  #   logwt <- rep(NA, nsamp)
  #   theta.tmp <- mvrnorm(nsamp, mu.theta, solve(Q.theta))
  #   for(i in 1:nsamp){
  #     logwt[i] <- F.logwt(theta.tmp[i,], spde, mu.theta, Q.theta, M)
  #   }
  #
  #   #weights to apply to each posterior sample of theta
  #   wt.tmp <- exp(logwt - max(logwt))
  #   wt <- wt.tmp/(sum(wt.tmp))
  #
  #
  #   ## Create index vectors
  #   K <- (dim(theta.tmp)[2] - 1)/2
  #   n.mesh <- mesh$n
  #   ind_beta <- list()
  #   for(k in 1:K){
  #     ind_beta[[k]] <- 1:n.mesh + (k-1)*n.mesh
  #   }
  #
  #   ## Compute cross-products for single session
  #   A.lst <- vector("list", K)
  #   for(k in 1:K){
  #     A.lst[[k]] <- A
  #   }
  #   Amat.tot <- bdiag(A.lst)
  #
  #   Xcros.all <- Xycros.all <- vector("list", M)
  #   for(m in 1:M){
  #     y_vec <- result[[m]]$y
  #     X_list <- result[[m]]$X
  #
  #     Xmat <- X_list[[1]]%*%Amat.tot
  #     Xcros.all[[m]] <- crossprod(Xmat)
  #     Xycros.all[[m]] <- crossprod(Xmat, y_vec)
  #   }
  #
  #   #get posterior quantities of beta, conditional on a value of theta
  #   if(is.null(no_cores)) {
  #     parallel <- FALSE
  #     beta.posteriors <- apply(theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, Xcros.all, Xycros.all, contrasts=contrasts, thresholds=thresholds, thresholds.contr=thresholds.contr, type=type, type.contr=type.contr, alpha=alpha, alpha.contr=alpha.contr, ind_beta=ind_beta)
  #   } else {
  #     ## Not sure how to change the following
  #     max_no_cores <- min(detectCores() - 1, 25)
  #     no_cores <- min(max_no_cores, no_cores)
  #     cl <- makeCluster(no_cores)
  #     beta.posteriors <- parApply(cl, theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, K=K, M, Xcros.all, Xycros.all, thresholds=thresholds, alpha=alpha, ind_beta=ind_beta)
  #     stopCluster(cl)
  #   }
  #
  #
  # ### Computing posterior quantities of beta averaged over subjects (summing over theta)
  # #organize samples
  # U <- length(thresholds)
  # mu.tot <- matrix(nrow=K*n.mesh, ncol=nsamp)
  # F.tot <- rep(list(rep(list(matrix(nrow=n.mesh, ncol=nsamp)), K)), U) #for each activation threshold and task, a Vx50 matrix
  # for(itheta in 1:nsamp){
  #   mu.tot[,itheta] <- beta.posteriors[[itheta]]$mean$mu
  #   for(u in 1:U){
  #     for(k in 1:K){
  #       F.tot[[u]][[k]][,itheta] <- beta.posteriors[[itheta]]$mean$F[[u]][,k]
  #     }
  #   }
  # }
  # ## Sum over samples using weights
  # betas.all <- matrix(0, nrow=n.mesh, ncol=K)
  # probs.all <- active.all <- array(0, dim=c(n.mesh, K, U)) #last dimension is for different activation thresholds
  #
  # #posterior mean
  # beta.pop <- as.vector(mu.tot%*%wt)
  # for(k in 1:K){
  #     beta.pop.k <- beta.pop[ind_beta[[k]]]
  #     betas.all[,k] <- as.vector(beta.pop.k)
  # }
  # #posterior probabilities and excursion sets
  # for(u in 1:U){
  #   for(k in 1:K){
  #     F.pop.uk <- as.vector(F.tot[[u]][[k]]%*%wt)
  #     E.pop.uk <- rep(0, length(F.pop.uk))
  #     E.pop.uk[F.pop.uk > 1-alpha] <- 1
  #     probs.all[,k,u] <- as.vector(F.pop.uk)
  #     active.all[,k,u] <- as.vector(E.pop.uk)
  #   }
  # }
  #
  #
  #
  # ### Computing posterior quantities of contrasts (summing over theta)
  # if(is.null(contrasts) == FALSE){
  # #organize samples
  # betas.all.contr <- probs.all.contr <- active.all.contr <- list()
  #  for(i.contr in 1:length(contrasts)){
  #   K.contr <- dim(contrasts[[i.contr]])[1]/n.mesh
  #   U.contr <- length(thresholds.contr[[i.contr]])
  #   mu.tot.contr <- matrix(nrow=K.contr*n.mesh, ncol=nsamp)
  #   F.tot.contr <- rep(list(rep(list(matrix(nrow=n.mesh, ncol=nsamp)), K.contr)), U.contr) #for each activation threshold and contrast, a Vx50 matrix
  #
  #   for(itheta in 1:nsamp){
  #     mu.tot.contr[,itheta] <- beta.posteriors[[itheta]]$contr$mu[[i.contr]]
  #     for(u in 1:U.contr){
  #       for(k in 1:K.contr){
  #         F.tot.contr[[u]][[k]][,itheta] <- beta.posteriors[[itheta]]$contr$F[[i.contr]][[u]][,k]
  #       }
  #     }
  #   }
  #   ## Sum over samples using weights
  #   betas.all.c <- matrix(0, nrow=n.mesh, ncol=K.contr)
  #   probs.all.c <- active.all.c <- array(0, dim=c(n.mesh, K.contr, U.contr)) #last dimension is for different activation thresholds
  #
  #   #posterior mean
  #   beta.pop.contr <- as.vector(mu.tot.contr%*%wt)
  #   for(k in 1:K.contr){
  #     beta.pop.k <- beta.pop.contr[ind_beta[[k]]]
  #     betas.all.c[,k] <- as.vector(beta.pop.k)
  #   }
  #   #posterior probabilities and excursion sets
  #   for(u in 1:U.contr){
  #     for(k in 1:K.contr){
  #       F.pop.uk <- as.vector(F.tot.contr[[u]][[k]]%*%wt)
  #       E.pop.uk <- rep(0, length(F.pop.uk))
  #       E.pop.uk[F.pop.uk > 1 - alpha.contr[i.contr]] <- 1
  #       probs.all.c[,k,u] <- as.vector(F.pop.uk)
  #       active.all.c[,k,u] <- as.vector(E.pop.uk)
  #     }
  #   }
  #   betas.all.contr[[i.contr]] <- betas.all.c
  #   probs.all.contr[[i.contr]] <- probs.all.c
  #   active.all.contr[[i.contr]] <- active.all.c
  #  }
  # } else{
  #   betas.all.contr <- probs.all.contr <- active.all.contr <- NULL
  # }
  #
  # ### Save all results
  # result <- list(mean = list(beta_estimates = betas.all, ppm = probs.all, active = active.all), contrasts = list(beta_estimates = betas.all.contr, ppm = probs.all.contr, active = active.all.contr))
  #
  # class(result) <- "BayesGLM_group"
  #
  # return(result)
  #
  # }
  #
  #

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
