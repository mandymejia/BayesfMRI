# Computes Bayesian GLM posterior means, quantiles and excursions sets (areas of activation) for summaries or contrasts across sessions within a single multi-session model
#
# @param result An object of class BayesGLM containing the results of a multi-session analysis. Must contain INLA_result!  (See `help(BayesGLM)`.)
# @param contrasts (Optional) A list of vectors, each length J*K, specifying the contrast(s) of interest across sessions, where J is the number of sessions analyzed and K is the number of tasks.  See Details for more information. Default is to compute the average for each task across sessions.
# @param quantiles (Optional) Vector of posterior quantiles to return in addition to the posterior mean
# @param excursion_type (For inference only) The type of excursion function for the contrast (">", "<", "!="), or a vector thereof (each element corresponding to one contrast).  If none provided, no inference performed.
# @param gamma (For inference only) Activation threshold for the excursion set, or a vector thereof (each element corresponding to one contrast).
# @param alpha (For inference only) Significance level for activation for the excursion set
# @param no_cores (For inference only) The number of cores to use for sampling in parallel. If NULL, do not run in parallel.
# @param nsamp_theta Number of theta values to sample from posterior. Default is 50.
# @param nsamp_beta Number of beta vectors to sample conditional on each theta value sampled. Default is 100.
#
# @details Performs posterior inference on summaries or contrasts across the sessions in a multi-session Bayesian GLM.  Let J be the number of sessions analyzed in the
# multi-session Bayesian GLM, and let K be the number of tasks.  Computes the posterior mean of each contrast and produces excursions sets representing areas of
# activation (above a specified activation threshold gamma) for each contrast map based on the joint (across locations) posterior distribution.
#
# Each contrast vector specifies an across-session summary of interest.  For example, the contrast vector `rep(c(1/J, rep(0, J-1)), K)` represents the across-session average
# for the first task; `c( rep(c(1/J1, rep(0, J1-1)), K), rep(c(-1/J2, rep(0, J2-1)), K))` represents the difference between the first J1 sessions and the remaining J2 sessions (J1+J2=J) for the first task;
#`rep( c(1/J,-1/J,rep(0, K-2)) ,J)` represents the difference between the first two tasks, averaged over all sessions.
#
# @return A list containing...
# @export
# @importFrom INLA inla.spde2.matern inla.pardiso.check inla.setOption inla.qsample
# @importFrom MASS mvrnorm
# @importFrom Matrix bdiag crossprod
# @importFrom abind abind
# @import parallel
#
# @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
# summarize_sessions.cifti <- function(result,
#                                      contrasts = NULL,
#                                      quantiles = NULL,
#                                      excursion_type,
#                                      gamma = 0,
#                                      alpha = 0.05,
#                                      #no_cores = NULL,
#                                      nsamp_theta = 50,
#                                      nsamp_beta = 100){
#
#   flag <- inla.pardiso.check()
#   if(grepl('FAILURE',flag)) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO is required for computational efficiency. See inla.pardiso().')
#   inla.setOption(smtp='pardiso')
#
#   # Find the numnber of sessions and tasks
#   session_names <- names(result$design)
#   J <- length(session_names)
#   K <- ncol(result$design[[1]])
#
#   #### SET UP OR CHECK CONTRAST VECTOR(S)
#
#   # If contrasts is null, by default set up a contrast vector that will compute the average across subjects/sessions for each task
#   if(!is.null(contrasts) & !is.list(contrasts)) contrasts <- list(contrasts)
#   if(is.null(contrasts)) {
#     contrasts <- vector('list', length=K)
#     for(k in 1:K){
#       contrast_onesubj_k <- c(rep(0, k-1), 1/J, rep(0, K-k)) #(1/J, 0, ... 0) for k=1, (0, 1/J, 0, ..., 0) for k=2, ..., (0, ..., 0, 1/J) for k=K
#       contrast_allsubj_k <- rep(contrast_onesubj_k, J)
#       contrasts[[k]] <- contrast_allsubj_k
#     }
#   }
#
#   #Check that each contrast vector is numeric and length J*K
#   num_contrasts <- length(contrasts)
#   length_each_contrast <- sapply(contrasts, length)
#   class_each_contrast <- sapply(contrasts, is.numeric)
#   if(any(length_each_contrast != J*K)) stop('each contrast vector must be of length J*K')
#   if(any(!class_each_contrast)) stop('each contrast vector must be numeric, but at least one is not')
#
#   #Check excursions settings
#   if(exists(excursion_type)) {
#     do_excur <- TRUE
#     if(length(excursion_type) == 1) excursion_type <- rep(excursion_type, num_contrasts)
#     if(length(gamma) == 1) gamma <- rep(gamma, num_contrasts)
#     if(length(alpha) != 1) stop('alpha must be a scalar value')
#     if(length(gamma) != num_contrasts) stop('Length of gamma must match number of contrasts or be equal to one.')
#     if(length(excursion_type) != num_contrasts) stop('Length of excursion_type must match number of contrasts or be equal to one.')
#   } else {
#     excursion_type <- 'none'
#   }
#
#   #Check quantiles argument
#   if(!is.null(quantiles)){
#     if(any(quantiles > 1 | quantiles < 0)) stop('All elements of quantiles argument must be between 0 and 1.')
#   }
#
#   #### LOOP OVER HEMISPHERES
#
#   do_left <- do_right <- TRUE
#   if(is.null(result$GLMs_Bayesian$cortex_left)) do_left <- FALSE
#   if(is.null(result$GLMs_Bayesian$cortex_right)) do_right <- FALSE
#
#   #put each hemisphere results into a list to create cifti objects later
#   posterior_means <- vector('list', length=2)
#   posterior_quantiles <- vector('list', length=2)
#   posterior_probs <- vector('list', length=2)
#   posterior_active <- vector('list', length=2)
#
#   for(h in 1:2){
#
#     ### CALL summarize_sessions()
#
#   } #end loop over hemispheres
#
#
#   ## Make cifti objects
#
#   posterior_means_cifti <- cifti_make(cortex_left = posterior_means[[1]],
#                                       cortex_right = posterior_means[[2]])
#                                       #surf_left = result_svh$betas_Bayesian$LR$SURF_LEFT$surface)
#   if(do_excur){
#     posterior_probs_cifti <- cifti_make(cortex_left = posterior_probs[[1]],
#                                         cortex_right = posterior_probs[[2]])
#     posterior_active_cifti <- cifti_make(cortex_left = posterior_active[[1]],
#                                          cortex_right = posterior_active[[2]])
#                                          #surf_left = result_svh$betas_Bayesian$LR$SURF_LEFT$surface)
#   }
#   if(num_quantiles > 0){
#     posterior_quantiles_cifti <- vector('list', num_quantiles)
#     names(posterior_quantiles_cifti) <- paste0('q',quantiles)
#     for(iq in 1:num_quantiles){
#       posterior_quantiles_cifti[[1]] <- cifti_make(cortex_left = posterior_quantiles[[1]][[iq]],
#                                                    cortex_right = posterior_quantiles[[2]][[iq]])
#     }
#   }
#
#
#   ### Save all results
#   result <- list(mean = posterior_means_cifti,
#                  quantiles = posterior_quantiles_cifti,
#                  PPM = posterior_probs_cifti,
#                  active = posterior_active_cifti)
#
#   class(result) <- "BayesGLM_group"
#
#   return(result)
#
# }

#' Summarize sessions
#' 
#' Computes Bayesian GLM posterior means, quantiles and excursions sets 
#'  (areas of activation) for summaries or contrasts across sessions within a 
#'  single multi-session model
#' 
#' Performs posterior inference on summaries or contrasts across the sessions 
#'  in a multi-session Bayesian GLM.  Let J be the number of sessions analyzed 
#'  in the multi-session Bayesian GLM, and let K be the number of tasks. 
#'  Computes the posterior mean of each contrast and produces excursions sets 
#'  representing areas of activation (above a specified activation threshold 
#'  gamma) for each contrast map based on the joint (across locations) posterior
#'  distribution.
#'
#' Each contrast vector specifies an across-session summary of interest. 
#'  The default is one contrast vector for each task representing the across-session average.
#'  Examples of contrast vectors include:
#'  \code{rep(c(1/J, rep(0, J-1)), K)} represents the across-session average for the first task.
#'  \code{rep( c(1/J,-1/J,rep(0, K-2)) ,J)} represents the difference between the first two tasks, averaged over all sessions.
#' 
#' @inheritSection INLA_Description INLA Requirement
#' 
#' @param result An object of class BayesGLM containing the results of a 
#'  multi-session analysis. Must contain INLA_result!  (See `help(BayesGLM)`.)
#' @param contrasts (Optional) A list of vectors, each length J*K, specifying 
#'  the contrast(s) of interest across sessions, where J is the number of 
#'  sessions analyzed and K is the number of tasks.  See Details for more 
#'  information. Default is to compute the average for each task across sessions.
#' @param quantiles (Optional) Vector of posterior quantiles to return in 
#'  addition to the posterior mean
#' @param excursion_type (For inference only) The type of excursion function 
#'  for the contrast (">", "<", "!="), or a vector thereof (each element 
#'  corresponding to one contrast).  If none provided, no inference performed.
#' @param gamma (For inference only) Activation threshold for the excursion set,
#'  or a vector thereof (each element corresponding to one contrast).
#' @param alpha (For inference only) Significance level for activation for the 
#'  excursion set
# @param no_cores (For inference only) The number of cores to use for sampling in parallel. If NULL, do not run in parallel.
#' @param nsamp_theta Number of theta values to sample from posterior. Default 
#'  is 50.
#' @param nsamp_beta Number of beta vectors to sample conditional on each theta
#'  value sampled. Default is 100.
#'
#' @return A list containing...
#' 
#' @importFrom INLA inla.spde2.matern inla.pardiso.check inla.setOption inla.qsample
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag crossprod
#' @importFrom abind abind
#' @import parallel
#' 
#' @export
summarize_sessions <- function(result,
                               contrasts = NULL,
                               quantiles = NULL,
                               excursion_type,
                               gamma = 0,
                               alpha = 0.05,
                               #no_cores = NULL,
                               nsamp_theta = 50,
                               nsamp_beta = 100){

  flag <- inla.pardiso.check()
  if(grepl('FAILURE',flag)) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO is required for computational efficiency. See inla.pardiso().')
  inla.setOption(smtp='pardiso')

  # Find the numnber of sessions and tasks
  session_names <- result$session_names # TO DO: have BayesGLM_surface return the design matrix.
  J <- length(session_names)
  K <- length(result$beta_names)

  #### SET UP OR CHECK CONTRAST VECTOR(S)

  # If contrasts is null, by default set up a contrast vector that will compute the average across subjects/sessions for each task
  if(!is.null(contrasts) & !is.list(contrasts)) contrasts <- list(contrasts)
  if(is.null(contrasts)) {
    contrasts <- vector('list', length=K)
    for(k in 1:K){
      contrast_onesubj_k <- c(rep(0, k-1), 1/J, rep(0, K-k)) #(1/J, 0, ... 0) for k=1, (0, 1/J, 0, ..., 0) for k=2, ..., (0, ..., 0, 1/J) for k=K
      contrast_allsubj_k <- rep(contrast_onesubj_k, J)
      contrasts[[k]] <- contrast_allsubj_k
    }
  }

  #Check that each contrast vector is numeric and length J*K
  num_contrasts <- length(contrasts)
  length_each_contrast <- sapply(contrasts, length)
  class_each_contrast <- sapply(contrasts, is.numeric)
  if(any(length_each_contrast != J*K)) stop('each contrast vector must be of length J*K')
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


  # Mesh and SPDE object
  mesh <- result$mesh
  spde <- inla.spde2.matern(mesh)
  mask <- result$mask

  #### GET THETA POSTERIOR FROM MULTI-SESSION MODEL

  Q.theta <- result$Q.theta
  mu.theta <- result$mu.theta

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

  Xcros <- Xycros <- vector('list', J)
  y_vec <- result$y
  X_list <- result$X #%*%Amat.tot #for multi-session data, make X a block diagonal matrix
  for(mm in 1:J){
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
        beta.posteriors[[isamp]] <- beta.posterior.thetasamp(theta.samp[,isamp],
                                                             spde=spde,
                                                             Xcros=Xcros,
                                                             Xycros=Xycros,
                                                             contrasts=contrasts,
                                                             quantiles=quantiles,
                                                             excursion_type=excursion_type,
                                                             gamma=gamma,
                                                             alpha=alpha,
                                                             nsamp_beta=nsamp_beta) #20 sec/iter
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

  return(list(
    posterior_means = mu.all,
    posterior_quantiles = quantiles.all,
    posterior_probs = F.all,
    posterior_active = active.all
  ))

}

