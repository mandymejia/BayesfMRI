#' Group-level Bayesian GLM
#'
#' Performs group-level Bayesian GLM estimation and inference using the joint
#'  approach described in Mejia et al. (2020)
#'
#' Each contrast vector specifies a group-level summary of interest. Let M be
#'  the number of subjects and K be the number of tasks. For example, the
#'  contrast vector `rep(rep(c(1/(M*num_sessions),rep(0, K-1)),num_sessions),M)` represents the group average
#'  for the first task for M subjects; `c(rep(rep(c(1/(M1*num_sessions),rep(0, K-1)),num_sessions),M1),rep(rep(c(-1/(M2*num_sessions),rep(0, K-1)),num_sessions),M2))`
#'  represents the difference between the first M1 subjects and the remaining M2
#'  subjects (M1+M2=M) for the first task; `rep(rep(c(1/(M*num_sessions),-1/(M*num_sessions),rep(0, K-2)), num_sessions),M)`
#'  represents the difference between the first two tasks, averaged over all
#'  subjects.
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param results Either (1) a list of length M of objects of class BayesGLM,
#'  or (2) a character vector of length M of file names output from the BayesGLM function.
#'  M is the number of subjects.
#' @param contrasts (Optional) A list of vectors, each length `M * K * num_sessions`, specifying the contrast(s)
#'  of interest across subjects, where M is the number of subjects and K is the number of tasks.
#'  See Details for more information. Default is to compute the average for each task across subjects.
#' @param quantiles (Optional) Vector of posterior quantiles to return in addition to the posterior mean
#' @param excursion_type (For inference only) The type of excursion function for the contrast (">", "<", "!="),
#'  or a vector thereof (each element corresponding to one contrast).  If NULL, no inference performed.
#' @param gamma (For inference only) List of vectors of activation thresholds for the excursion set (each element corresponding to one contrast). Remember that if a contrast is not specified, the average is found.
#' @param alpha (For inference only) Significance level for activation for the excursion set, or a vector thereof (each element corresponding to one contrast).
#' @param nsamp_theta Number of theta values to sample from posterior. Default is 50.
#' @param nsamp_beta Number of beta vectors to sample conditional on each theta value sampled. Default is 100.
#' @param no_cores The number of cores to use for sampling betas in parallel. If NULL, do not run in parallel.
#' @inheritParams verbose_Param_direct_TRUE
#' @param use_EM (logical) Should the method for the EM algorithm be used?
#'
#' @return A list containing the estimates, PPMs and areas of activation for each contrast.
#'
#' @importFrom INLA inla.spde2.matern inla.spde.make.A
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag crossprod
#'
#' @export
BayesGLM2 <- function(results,
                           contrasts = NULL,
                           quantiles = NULL,
                           excursion_type=NULL,
                           gamma = list(0),
                           alpha = 0.05,
                           nsamp_theta = 50,
                           nsamp_beta = 100,
                           no_cores = NULL,
                           verbose = TRUE,
                           use_EM = FALSE){

  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("`BayesGLM2` requires the `abind` package. Please install it.", call. = FALSE)
  }

  # Check to see that the INLA package is installed
  check_BayesGLM(require_PARDISO=TRUE)

  #Check if results are model objects or file paths

  result_class <- sapply(results, class)
  names(result_class) <- NULL
  problem <- FALSE
  is_RDS <- function(x){
    if(!is.character(x)) {
      return(FALSE)
    } else {
      # Damon: Replaced this line to reduce number of package dependencies:
      #   stringr::str_sub(x, -3, -1)
      # https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
      last3 <- sub(".*(?=.{3}$)", "", x, perl=TRUE)
      return(last3 %in% c('rds','RDS','Rds'))
    }
  }
  if(result_class[1] == 'character'){
    if(any(!is_RDS(results))) problem <- TRUE
  } else if(result_class[1] == 'BayesGLM'){
    if(any(result_class != 'BayesGLM')) problem <- TRUE
  } else {
    problem <- TRUE
  }
  if(problem) stop('All elements of results argument must be BayesGLM or a character file path to an RDS object.')

  #if using files, grab the first model result for checks
  M <- length(results)
  use_files <- (result_class[1]=='character')
  use_objects <- (result_class[1]=='BayesGLM')
  if(use_files) {
    fnames <- results
    results <- vector('list', length=M)
    results[[1]] <- readRDS(fnames[1])
    if(class(results[[1]]) == "BayesGLM_cifti") {
      if(use_EM) {
        which_BayesGLM <- which(sapply(results[[1]]$GLMs_EM,class) == "BayesGLM")
        results[[1]] <- results[[1]]$GLMs_EM[[which_BayesGLM]]
      }
      if(!use_EM) {
        which_BayesGLM <- which(sapply(results[[1]]$GLMs_Bayesian,class) == "BayesGLM")
        results[[1]] <- results[[1]]$GLMs_Bayesian[[which_BayesGLM]]
      }
    }
    if(class(results[[1]]) != 'BayesGLM') stop("Each RDS file in results argument must contain an object of class BayesGLM")
  }


  #Check that subject-level models are single-session models
  num_sessions <- length(results[[1]]$session_names)
  # if(num_sessions>1) stop('This function is currently only applicable to results of single-session modeling at subject level.')
  # We actually can't really use the averages here because we don't have a
  # corresponding design matrix and response to calculate the posteriors for beta.
  # Instead, we can get this working for multi-session data. To begin, we will
  # assume the same number of sessions for each subject, but this can be changed
  # later.
  # if(num_sessions > 1) {
  #   message("Currently, multi-session analysis is only possible using the averages across sessions.")
  #   use_avg <- "matrix" %in% class(results[[1]]$avg_beta_estimates)
  #   if(!use_avg) stop("This function does not currently support group analysis on multiple sessions unless the analysis is done on subject-level averages across sessions.")
  #   num_sessions <- 1
  # }
  # use_avg <- FALSE

  #### SET UP OR CHECK CONTRAST VECTOR(S)

  #TO DO: Generalize to multi-session models (just need to expand the contrasts vector to summarize across sessions as well as subjects)

  # If contrasts is null, by default set up a contrast vector that will compute the average across subjects for each task
  K <- length(results[[1]]$beta_names) #number of tasks
  beta_names <- results[[1]]$beta_names
  if(!is.null(contrasts) & !is.list(contrasts)) contrasts <- list(contrasts)
  if(is.null(contrasts)) {
    if(verbose) cat('Setting up contrast vectors to compute the average across subjects for each task. If other contrasts are desired, provide a contrasts argument.\n')
    contrasts <- vector('list', length=K)
    names(contrasts) <- paste0(beta_names,'_avg')
    for(k in 1:K){
      # Old, broken, fixed?
      contrast_onesubj_k <- rep(c(rep(0, k-1), 1/(num_sessions*M), rep(0, K-k)),num_sessions) #(1/J, 0, ... 0) for k=1, (0, 1/J, 0, ..., 0) for k=2, ..., (0, ..., 0, 1/J) for k=K, for each session
      contrast_allsubj_k <- rep(contrast_onesubj_k, M)
      contrasts[[k]] <- contrast_allsubj_k
    }
  }

  #Check that each contrast vector is numeric and length J*K
  num_contrasts <- length(contrasts)
  length_each_contrast <- sapply(contrasts, length)
  class_each_contrast <- sapply(contrasts, is.numeric)
  if(any(length_each_contrast != M*K*num_sessions)) stop('each contrast vector must be of length M*K*num_sessions')
  if(any(!class_each_contrast)) stop('each contrast vector must be numeric, but at least one is not')

  #Check excursions settings
  if(is.null(excursion_type)) do_excur <- FALSE else do_excur <- TRUE
  if(do_excur) {
    if(length(excursion_type) == 1) excursion_type <- rep(excursion_type, num_contrasts)
    if(length(gamma) == 1) gamma <- rep(gamma, num_contrasts)
    if(length(alpha) == 1) alpha <- rep(alpha, num_contrasts)
    if(length(gamma) != num_contrasts) stop('Length of gamma must match number of contrasts or be equal to one.')
    if(length(alpha) != num_contrasts) stop('Length of alpha must match number of contrasts or be equal to one.')
    if(length(excursion_type) != num_contrasts) stop('Length of excursion_type must match number of contrasts or be equal to one.')
  } else {
    excursion_type <- 'none'
  }

  # Mesh and SPDE object
  mesh <- results[[1]]$mesh
  if("Amat" %in% names(mesh) & "spde" %in% names(mesh)){
    spde <- mesh$spde
    Amat <- mesh$Amat
  } else {
    spde <- inla.spde2.matern(mesh)
    Amat <- inla.spde.make.A(mesh) #Psi_{km} (for one task and subject, a VxN matrix, V=num_vox, N=num_mesh)
    Amat <- Amat[mesh$idx$loc,]
  }

  Amat.tot <- bdiag(rep(list(Amat), K)) #Psi_m from paper (VKxNK)

  #Check quantiles argument
  if(!is.null(quantiles)){
    if(any(quantiles > 1 | quantiles < 0)) stop('All elements of quantiles argument must be between 0 and 1.')
  }

  # Collecting theta posteriors from subject models
  if(!use_EM) {
    Qmu_theta <- Q_theta <- 0
  }
  if(use_EM) mu_theta <- 0

  # Collecting X and y cross-products from subject models (for posterior distribution of beta)
  Xcros.all <- Xycros.all <- vector("list", M)

  for(m in 1:M){

    if(use_files & (m > 1)){
      results[[m]] <- readRDS(fnames[m])
      if(class(results[[m]]) == "BayesGLM_cifti") {
        if(use_EM) {
          which_BayesGLM <- which(sapply(results[[m]]$GLMs_EM,class) == "BayesGLM")
          results[[m]] <- results[[m]]$GLMs_EM[[which_BayesGLM]]
        }
        if(!use_EM) {
          which_BayesGLM <- which(sapply(results[[m]]$GLMs_Bayesian,class) == "BayesGLM")
          results[[m]] <- results[[m]]$GLMs_Bayesian[[which_BayesGLM]]
        }
      }
      if(class(results[[m]]) != 'BayesGLM') stop("Each RDS file in results argument must contain an object of class BayesGLM")
      # use_avg_m <- "matrix" %in% class(results[[m]]$avg_beta_estimates)
      num_sessions_m <- length(results[[m]]$session_names)
      # if(use_avg & use_avg_m) num_sessions_m <- 1
      if(num_sessions_m != num_sessions) stop(paste("Group modeling currently only supports an equal number of sessions across all subjects. Subject 1 has", num_sessions, "sessions, but subject", m, "has",paste0(length(results[[m]]$session_names),".")))
    }
    if(m == 1) num_sessions_m <- num_sessions
    # use_avg_m <- "matrix" %in% class(results[[m]]$avg_beta_estimates)

    #Check match of beta names
    beta_names_m <- results[[m]]$beta_names
    if(!all.equal(beta_names_m, beta_names, check.attributes=FALSE)) stop('All subjects must have the same tasks in the same order.')

    #Check that mesh has same neighborhood structure
    faces_m <- results[[m]]$mesh$faces
    if(!all.equal(faces_m, mesh$faces, check.attribute=FALSE)) stop(paste0('Subject ',m,' does not have the same mesh neighborhood structure as subject 1. Check meshes for discrepancies.'))

    #Check that model is single session or average
    # if(num_sessions_m>1 & !use_avg_m)
    # stop('This function is currently only applicable to results of single-session or average modeling at subject level.')

    #Collect posterior mean and precision of hyperparameters
    mu_theta_m <- results[[m]]$mu.theta
    if(!use_EM) {
      Q_theta_m <- results[[m]]$Q.theta
      #iteratively compute Q_theta and mu_theta (mean and precision of q(theta|y))
      Qmu_theta <- Qmu_theta + as.vector(Q_theta_m%*%mu_theta_m)
      Q_theta <- Q_theta + Q_theta_m
      rm(mu_theta_m, Q_theta_m)
    }
    if(use_EM) {
      mu_theta <- mu_theta + mu_theta_m
      rm(mu_theta_m)
    }

    #compute Xcros = Psi'X'XPsi and Xycros = Psi'X'y (all these matrices for a specific subject m)
    y_vec <- results[[m]]$y
    X_list <- results[[m]]$X
    if(length(X_list) > 1) {
      n_sess <- length(X_list)
      X_list <- Matrix::bdiag(X_list)
      Amat.final <- Matrix::bdiag(rep(list(Amat.tot),n_sess))
    } else {
      X_list <- X_list[[1]]
      Amat.final <- Amat.tot
    }
    Xmat <- X_list%*%Amat.final
    Xcros.all[[m]] <- Matrix::crossprod(Xmat)
    Xycros.all[[m]] <- Matrix::crossprod(Xmat, y_vec)

    if(m > 1) results[[m]] <- c(0) #to save memory
  }
  if(!use_EM) {
    mu_theta <- solve(Q_theta, Qmu_theta) #mu_theta = posterior mean of q(theta|y) (Normal approximation) from paper, Q_theta = posterior precision

    #### DRAW SAMPLES FROM q(theta|y)

    #theta.tmp <- mvrnorm(nsamp_theta, mu.theta, solve(Q.theta))
    if(verbose) cat(paste0('Sampling ',nsamp_theta,' posterior samples of thetas \n'))
    theta.samp <- inla.qsample(n=nsamp_theta, Q = Q_theta, mu = mu_theta)

    #### COMPUTE WEIGHT OF EACH SAMPLES FROM q(theta|y) BASED ON PRIOR

    if(verbose) cat('Computing weights for each theta sample \n')
    logwt <- rep(NA, nsamp_theta)
    for(i in 1:nsamp_theta){ logwt[i] <- F.logwt(theta.samp[,i], spde, mu_theta, Q_theta, M) }
    #weights to apply to each posterior sample of theta
    wt.tmp <- exp(logwt - max(logwt))
    wt <- wt.tmp/(sum(wt.tmp))
  }
  if(use_EM) {
    mu_theta <- mu_theta / M
    theta.samp <- as.matrix(mu_theta)
    wt <- 1
  }

  #get posterior quantities of beta, conditional on a value of theta
  if(verbose) cat(paste0('Sampling ',nsamp_beta,' betas for each value of theta \n'))
  if(is.null(no_cores)){
    #6 minutes in simulation
    beta.posteriors <- apply(theta.samp,
                             MARGIN=2,
                             FUN=beta.posterior.thetasamp,
                             spde=spde,
                             Xcros = Xcros.all,
                             Xycros = Xycros.all,
                             contrasts=contrasts,
                             quantiles=quantiles,
                             excursion_type=excursion_type,
                             gamma=gamma,
                             alpha=alpha,
                             nsamp_beta=nsamp_beta)
  } else {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("`BayesGLM2` requires the `parallel` package. Please install it.", call. = FALSE)
    }

    #2 minutes in simulation (4 cores)
    max_no_cores <- min(parallel::detectCores() - 1, 25)
    no_cores <- min(max_no_cores, no_cores)
    cl <- parallel::makeCluster(no_cores)
    beta.posteriors <- parallel::parApply(cl, theta.samp,
                                          MARGIN=2,
                                          FUN=beta.posterior.thetasamp,
                                          spde=spde,
                                          Xcros = Xcros.all,
                                          Xycros = Xycros.all,
                                          contrasts=contrasts,
                                          quantiles=quantiles,
                                          excursion_type=excursion_type,
                                          gamma=gamma,
                                          alpha=alpha,
                                          nsamp_beta=nsamp_beta)
    parallel::stopCluster(cl)
  }

  ## Sum over samples using weights

  if(verbose) cat('Computing weighted summaries over beta samples \n')

  ## Posterior mean of each contrast
  betas.all <- lapply(beta.posteriors, function(x) return(x$mu))
  betas.wt <- mapply(function(x, a){return(x*a)}, betas.all, wt, SIMPLIFY=FALSE) #apply weight to each element of betas.all (one for each theta sample)
  betas.summ <- apply(abind::abind(betas.wt, along=3), MARGIN = c(1,2), sum)  #N x L (# of contrasts)

  ## Posterior quantiles of each contrast
  num_quantiles <- length(quantiles)
  if(num_quantiles > 0){
    quantiles.summ <- vector('list', num_quantiles)
    names(quantiles.summ) <- quantiles
    for(iq in 1:num_quantiles){
      quantiles.all_iq <- lapply(beta.posteriors, function(x) return(x$quantiles[[iq]]))
      betas.wt_iq <- mapply(function(x, a){return(x*a)}, quantiles.all_iq, wt, SIMPLIFY=FALSE) #apply weight to each element of quantiles.all_iq (one for each theta sample)
      quantiles.summ[[iq]] <- apply(abind::abind(betas.wt_iq, along=3), MARGIN = c(1,2), sum)  #N x L (# of contrasts)
    }
  } else {
    quantiles.summ <- NULL
  }

  ## Posterior probabilities and activations
  if(do_excur){
    ppm.all <- lapply(beta.posteriors, function(x) return(x$F))
    ppm.wt <- mapply(function(x, a){sapply(x,function(xx,aa) xx*aa, aa = a, simplify = F)}, ppm.all, wt, SIMPLIFY=FALSE) #apply weight to each element of ppm.all (one for each theta sample)
    # ppm.wt <- mapply(function(x, a){return(x*a)}, ppm.all, wt, SIMPLIFY=FALSE) #apply weight to each element of ppm.all (one for each theta sample)
    ppm.summ <- Reduce(function(x,y) {
      output_list <- vector("list",num_contrasts)
      for(l in 1:num_contrasts) {
        output_list[[l]] <- x[[l]] + y[[l]]
      }
      return(output_list)
    }, ppm.wt)
    # ppm.summ <- apply(abind::abind(ppm.wt, along=3), MARGIN = c(1,2), sum) #N x L (# of contrasts)
    active <- mapply(function(ppm,al){
      output_matrix <- array(0,dim = dim(as.matrix(ppm)))
      output_matrix[c(ppm) > (1 - al)] <- 1
      return(output_matrix)
    }, ppm = ppm.summ, al = alpha, SIMPLIFY = F)
    # active <- array(0, dim=dim(ppm.summ))
    # for(l in 1:num_contrasts){
    #   active[ppm.summ[,l] > (1-alpha[l]),l] <- 1
    # }
  } else {
    ppm.summ <- active <- NULL
  }

  ### Save all results
  result <- list(estimates = betas.summ,
                 quantiles = quantiles.summ,
                 ppm = ppm.summ,
                 active = active,
                 contrasts = contrasts,
                 gamma = gamma,
                 alpha = alpha,
                 nsamp_theta = nsamp_theta,
                 nsamp_beta = nsamp_beta,
                 Amat = Amat)

  return(result)

}


#' Beta posterior theta sampling
#'
#' Internal function used in joint approach to group-analysis
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param theta A single sample of theta (hyperparameters) from q(theta|y)
#' @param spde A SPDE object from inla.spde2.matern() function.
#' @param Xcros A crossproduct of design matrix.
#' @param Xycros A crossproduct of design matrix and BOLD y.
#' @param contrasts A list of vectors of length M*K specifying the contrasts of interest.
#' @param quantiles Vector of posterior quantiles to return in addition to the posterior mean
#' @param excursion_type Vector of excursion function type (">", "<", "!=") for each contrast
#' @param gamma List of vector of activation thresholds. Each list element should correspond to a contrast
#' @param alpha Significance level for activation for the excursion sets
#' @param nsamp_beta The number of samples to draw from full conditional of beta given the current value of theta (p(beta|theta,y))
#'
#' @importFrom excursions excursions.mc
#' @importFrom Matrix Diagonal
#' @importFrom INLA inla.spde2.precision inla.qsample inla.qsolve
#'
#' @return A list containing...
#'
#' @keywords internal
beta.posterior.thetasamp <- function(theta,
                                     spde,
                                     Xcros,
                                     Xycros,
                                     contrasts,
                                     quantiles,
                                     excursion_type,
                                     gamma,
                                     alpha,
                                     nsamp_beta=100){

  n.mesh <- spde$n.spde

  # print('contructing joint precision')
  prec.error <- exp(theta[1])
  theta_spde <- matrix(theta[-1], nrow=2) #2xK matrix of the hyperparameters (2 per task)
  K <- ncol(theta_spde)
  M <- length(Xcros)

  use_EM <- all(sapply(c("M0","M1","M2"), function(x) x %in% names(spde)))

  #contruct prior precision matrix for beta, Q_theta for given sampled value of thetas
  # For EM
  if(use_EM) {
    Q.beta <- apply(theta_spde,2,function(theta_k) {
      theta_k <- exp(theta_k)^2
      out <- theta_k[1] * (theta_k[2]^2 * spde$M0 + theta_k[2] * spde$M1 + spde$M2)
      return(out)
    })
  }
  # For INLA
  if(!use_EM) {
    Q.beta <- list()
    for(k in 1:K) {
      theta_k <- theta_spde[,k] #theta[(2:3) + 2*(k-1)] #1:2, 2:3, 4:5, ...
      Q.beta[[k]] <- inla.spde2.precision(spde, theta = theta_k) # prior precision for a single task k
    }
  }
  Q <- Matrix::bdiag(Q.beta) #Q_theta in the paper

  N <- dim(Q.beta[[1]])[1] #number of mesh locations
  if(N != n.mesh) stop('Length of betas does not match number of vertices in mesh. Inform developer.')

  beta.samples <- NULL
  #~5 seconds per subject with PARDISO
  # print('Looping over subjects or sessions')
  num_sessions <- 1
  Q_mm <- Q
  for(mm in 1:M){
    if(nrow(Q) != nrow(Xcros[[mm]])) {
      num_sessions <- nrow(Xcros[[mm]]) / nrow(Q)
      Q_mm <- Matrix::bdiag(rep(list(Q),num_sessions))
    }
    #compute posterior mean and precision of beta|theta
    Q.m <- prec.error*Xcros[[mm]] + Q_mm #Q_m in paper
    mu.m <- inla.qsolve(Q.m, prec.error*Xycros[[mm]]) #NK x 1 -- 2 minutes, but only 2 sec with PARDISO!  #mu_m in paper

    #draw samples from pi(beta_m|theta,y)
    beta.samp.m <- inla.qsample(n = nsamp_beta, Q = Q.m, mu = mu.m) #NK x nsamp_beta  -- 2 minutes, but only 2 sec with PARDISO!

    #concatenate samples over models
    beta.samples <- rbind(beta.samples, beta.samp.m) #will be a (N*K*M*num_sessions) x nsamp_beta matrix

  }

  if(excursion_type[1] == 'none') do_excur <- FALSE else do_excur <- TRUE

  # Loop over contrasts
  num_contrasts <- length(contrasts)
  # mu.contr <- sapply(seq(num_contrasts), function(l) {
  #   matrix(NA, nrow = n.mesh, ncol = length(gamma[[l]]))
  # }, simplify = F)
  mu.contr <- matrix(NA, nrow=n.mesh, ncol=num_contrasts)
  if(do_excur) {
    F.contr <- sapply(seq(num_contrasts), function(l) {
      matrix(NA, nrow = n.mesh, ncol = length(gamma[[l]]))
    }, simplify = F)
  } else F.contr <- NULL
  if(!is.null(quantiles)){
    num_quantiles <- length(quantiles)
    quantiles.contr <- rep(list(mu.contr), num_quantiles)
    names(quantiles.contr) <- quantiles
  } else {
    num_quantiles <- 0
    quantiles.contr <- NULL
  }
  for(l in 1:num_contrasts){

    #Construct "A" matrix from paper (linear combinations)
    ctr.mat <- kronecker(t(contrasts[[l]]), Diagonal(n.mesh, 1))

    #beta.mean.pop.contr <- as.vector(ctr.mat%*%beta.mean.pop.mat)  # NKx1 or Nx1
    samples_l <- as.matrix(ctr.mat%*%beta.samples)  # N x nsamp_beta
    # For multiple sessions, add the beta samples for the sessions together
    # if(num_sessions > 1) {
    #   session_betas <- sapply(seq(num_sessions), function(sess) {
    #     session_inds <- seq(N) + (sess - 1)*N
    #     return(samples_l[session_inds,])
    #   }, simplify = FALSE)
    #   samples_l <- Reduce(`+`, session_betas)
    # }
    # mu.contr[[l]] <- rowMeans(samples_l) #compute mean over beta samples
    mu.contr[,l] <- rowMeans(samples_l) #compute mean over beta samples
    if(num_quantiles > 0){
      for(iq in 1:num_quantiles){
        quantiles.contr[[iq]][,l] <- apply(samples_l, 1, quantile, quantiles[iq])
      }
    }

    # Estimate excursions set for current contrast
    if(do_excur){
      excur_l <- sapply(gamma[[l]], function(each_gamma) {
        excur_lm <- excursions::excursions.mc(samples_l, u = each_gamma, type = excursion_type[l], alpha = alpha[l])
      }, simplify = FALSE)
      names(excur_l) <- as.character(gamma[[l]])
      F.contr[[l]] <- Reduce(cbind,sapply(excur_l,getElement,"F", simplify = F))
      # F.contr[,l] <- excur_l$F
    }
  }
  names(F.contr) <- names(contrasts)
  result <- list(mu = mu.contr, quantiles=quantiles.contr, F = F.contr)
  return(result)
}


#' F logwt
#'
#' Internal function used in joint approach to group-analysis for combining across models
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param theta A vector of hyperparameter values at which to compute the posterior log density
#' @param spde A SPDE object from inla.spde2.matern() function, determines prior precision matrix
#' @param mu_theta Posterior mean from combined subject-level models.
#' @param Q_theta Posterior precision matrix from combined subject-level models.
#' @param M Number of subjects
#' @return A list containing...
#'
#' @importFrom stats dgamma
#'
#' @keywords internal
F.logwt <- function(theta, spde, mu_theta, Q_theta, M){
  #mu_theta - posterior mean from combined subject-level models
  #Q_theta - posterior precision matrix from combined subject-level models
  #M - number of subjects
  a <- 1; b <- 5e-5
  n.spde <- (length(theta) - 1)/2
  mu.tmp <- spde$f$hyper$theta1$param[1:2]
  mu <- rep(mu.tmp, n.spde)
  Q.tmp <- matrix(spde$f$hyper$theta1$param[-(1:2)], 2, 2, byrow = TRUE)
  Q <- kronecker(diag(1, n.spde, n.spde), Q.tmp)

  ## Prior density
  pr.delta <- dgamma(exp(theta[1]), a, b, log = TRUE) #log prior density on residual precision
  pr.tk <- as.vector(-t(theta[-1] - mu)%*%Q%*%(theta[-1] - mu))/2 + log(det(Q))/2 - dim(Q)[1]*log(2*pi)/2 #joint log prior density on 2K spde parameters
  pr.theta <- pr.delta + pr.tk

  (1-M)*pr.theta
}

#' @rdname BayesGLM2
#' @export
BayesGLM_group <- function(
  results,
  contrasts = NULL,
  quantiles = NULL,
  excursion_type=NULL,
  gamma = 0,
  alpha = 0.05,
  nsamp_theta = 50,
  nsamp_beta = 100,
  no_cores = NULL,
  verbose = TRUE){

  BayesGLM2(
    results=results,
    contrasts=contrasts,
    quantiles=quantiles,
    excursion_type=excursion_type,
    gamma=gamma, alpha=alpha,
    nsamp_theta=nsamp_theta, nsamp_beta=nsamp_beta,
    no_cores=no_cores, verbose=verbose
  )
}

#' Group-level Bayesian GLM for subcortical data
#'
#' This is a wrapper function that applies group modeling to subcortical
#' results. This currently only works for output from the EM implementation of
#' the Bayesian GLM.
#'
#' @param results A character vector of length M of file names output from the
#'   BayesGLM_cifti function with subcortical results. M is the number of subjects.
#' @param contrasts (Optional) A list of vectors, each length `M * K * num_sessions`, specifying the contrast(s)
#'  of interest across subjects, where M is the number of subjects and K is the number of tasks.
#'  See Details for more information. Default is to compute the average for each task across subjects.
#' @param quantiles (Optional) Vector of posterior quantiles to return in addition to the posterior mean
#' @param excursion_type (For inference only) The type of excursion function for the contrast (">", "<", "!="),
#'  or a vector thereof (each element corresponding to one contrast).  If NULL, no inference performed.
#' @param gamma (For inference only) List of vectors of activation thresholds for the excursion set (each element corresponding to one contrast). Remember that if a contrast is not specified, the average is found.
#' @param alpha (For inference only) Significance level for activation for the excursion set, or a vector thereof (each element corresponding to one contrast).
#' @param nsamp_theta Number of theta values to sample from posterior. Default is 50.
#' @param nsamp_beta Number of beta vectors to sample conditional on each theta value sampled. Default is 100.
#' @inheritParams verbose_Param_direct_TRUE
#'
#' @return A list with length equal to the number of subcortical models run in each of the single-subject data cases
#' @export
BayesGLM2_vol <- function(results,
                          contrasts = NULL,
                          quantiles = NULL,
                          excursion_type=NULL,
                          gamma = list(0),
                          alpha = 0.05,
                          nsamp_theta = 50,
                          nsamp_beta = 100,
                          verbose = TRUE) {
  M <- length(results)
  first_result <- readRDS(results[1])
  which_BayesGLM <- which(sapply(first_result$GLMs_EM,class) == "BayesGLM")
  cifti_output <- first_result$betas_EM[[1]]
  # data_idx <- !is.na(cifti_output$data$subcort[,1])
  cifti_output$data$subcort[!is.na(cifti_output$data$subcort)] <- 0
  if(!3 %in% which_BayesGLM) stop("This function only works with subcortical results.")
  first_result <- first_result$GLMs_EM[[which_BayesGLM]]$EM_result_all
  num_regions <- length(first_result)
  for(reg in seq(num_regions)) {
    first_result[[reg]]$posterior_Sig_inv <- 0 #Done to save memory
  }
  results_out <- vector("list", length = M)
  names(results_out) <- names(first_result)

  for(reg in seq(num_regions)){
    cat("Evaluating group model for subcortical model",names(first_result)[reg],"\n")
    results_in <- vector("list",length = M)
    results_in[[1]] <- first_result[[reg]]
    class(results_in[[1]]) <- "BayesGLM"
    if(M > 1){
      for(m in 2:M) {
        results_in[[m]] <- readRDS(results[m])$GLMs_EM[[which_BayesGLM]]$EM_result_all[[reg]]
        results_in[[m]]$posterior_Sig_inv <- 0 # Done to save memory
        class(results_in[[m]]) <- "BayesGLM"
      }
    }
    results_out[[reg]] <- BayesGLM2(results = results_in,
                                    contrasts = contrasts,
                                    quantiles = quantiles,
                                    excursion_type=excursion_type,
                                    gamma = gamma,
                                    alpha = alpha,
                                    nsamp_theta = nsamp_theta,
                                    nsamp_beta = nsamp_beta,
                                    no_cores = NULL,
                                    verbose = verbose,
                                    use_EM = TRUE)
    results_out[[reg]]$cifti_estimate <- as.matrix(results_out[[reg]]$Amat %*% results_out[[reg]]$estimates)
    results_out[[reg]]$cifti_active <- sapply(results_out[[reg]]$active, function(x) as.matrix(results_out[[reg]]$Amat %*% x), simplify = F)
    results_out[[reg]]$cifti_idx <- unlist(first_result[[reg]]$mesh$idx)
    active_levels <- Reduce(cbind,sapply(results_out[[reg]]$cifti_active, function(x) apply(x,1,sum), simplify = F))
    active_levels[active_levels == 0] <- NA
    cifti_output$data$subcort[results_out[[reg]]$cifti_idx,] <- active_levels
  }
  return(results_out)
}
