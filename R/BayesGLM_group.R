#' Applies joint approach to group-level analysis to task fMRI data
#'
#' @param results Either (1) a list of length M of objects of class BayesGLM, or (2) a character vector of length M of file names output from the BayesGLM function. M is the number of subjects.
#' @param contrast A vector specifying the contrast of interest.  See Details for more information.
#' @param no_cores The number of cores to use for sampling in parallel
#'
#' @details The contrast vector specifies the group-level quantity of interest.  For example, the vector rep(1,M*K) would return the group average for each of K tasks;
#' the vector `c(rep(1,M1*K)`, `rep(-1,M2*K))` would return the difference between the average within two groups of size M1 and M2, respectively, for each of K tasks;
#' the vector `rep(rep(1,-1,0,0),each=V),M)` would return the difference between the first two tasks (of 4), averaged over all subjects.
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern
#'
#' @examples \dontrun{}
BayesGLM_group <- function(results, contrasts, no_cores=NULL){
  #check whether data is a list OR a session (for single-session analysis)
  #check whether each element of data is a session (use is.session)
  # V = number of data locations
  # T = length of time series for each session (vector)
  # K = number of tasks for each subject (each subject must have the same tasks)
  # M = number of subjects

  #check that each subject-level object is of class BayesGLM_obj
  #check that beta_names matches for all subjects
  #check that mesh (triangle) matches for all subjects

  #check that results argument is in one of the two correct formats
  #set up A matrix using the contrasts argument

  # #check that only mesh OR vertices+faces supplied
  # has_mesh <- !is.null(mesh)
  # has_verts_faces <- !is.null(vertices) & !is.null(faces)
  # has_howmany <- has_mesh + has_verts_faces
  # if(has_howmany != 1) stop('Must supply EITHER mesh OR vertices and faces.')
  #
  # #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  #
  # if(!is.list(data)) stop('I expect data to be a list, but it is not')
  #   data_classes <- sapply(data, 'class')
  # if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  # Find the numnber of subjects.
  subject_names <- names(data)
  M <- length(subject_names)

  V <- ncol(data[[1]]$BOLD) #number of data locations
  K <- ncol(data[[1]]$design) #number of tasks
  for(s in 1:M){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  }

  if(is.null(outfile)){
    warning('No value supplied for outfile, which is required for post-hoc group modeling.')
  }

  if(is.null(mesh)) {
    mesh <- make_mesh(vertices, faces)
  }


  spde <- inla.spde2.matern(mesh)

  # Collecting theta posteriors from subject models
  theta.sub <- NULL
  mu.theta.tmp <- Q.theta <- 0
  Xcros.all <- Xycros.all <- vector("list", M)

  #construct betas and repls objects
  replicates_list <- organize_replicates(n_sess=1, n_task=K, mesh=mesh)
  betas <- replicates_list$betas
  repls <- replicates_list$repls

  #organize the formula and data objects
  #formula <- make_formula(beta_names = names(betas), repl_names = names(repls), hyper_initial = c(-2,2))
  #formula <- as.formula(formula)

  beta_names <- names(betas)
  repl_names <- names(repls)
  n_beta <- length(names(betas))
  hyper_initial <- c(-2,2)
  hyper_initial <- rep(list(hyper_initial), n_beta)
  hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

  formula_vec <- paste0('f(',beta_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
  formula_vec <- c('y ~ -1', formula_vec)
  formula_str <- paste(formula_vec, collapse=' + ')
  formula <- as.formula(formula_str, env = globalenv())

  for(s in 1:M){

      #collect data and design matrices
      #extract and mask BOLD data for current session
      BOLD_s <- data[[s]]$BOLD

      #scale data to represent % signal change (or just center if scale=FALSE)
      BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale)
      design_s <- scale(data[[s]]$design, scale=FALSE) #center design matrix to eliminate baseline

      #regress nuisance parameters from BOLD data and design matrix
      if('nuisance' %in% names(data[[s]])){
        design_s <- data[[s]]$design
        nuisance_s <- data[[s]]$nuisance
        y_reg <- nuisance_regress(BOLD_s, nuisance_s)
        X_reg <- nuisance_regress(design_s, nuisance_s)
      } else {
        y_reg <- BOLD_s
        X_reg <- data[[s]]$design
      }

      #set up data and design matrix
      data_org <- organize_data(y_reg, X_reg)
      y_vec <- data_org$y
      X_list <- list(data_org$A)
      # names(X_list) <- session_names[s]

      # Computing cross-products for each subject
      Xcros.all[[mm]] <- crossprod(X_list[[1]])
      Xycros.all[[mm]] <- crossprod(X_list[[1]], y_vec)


      model_data <- make_data_list(y=y_vec, X=X_list, betas=betas, repls=repls)

      #estimate model using INLA
      INLA_result <- estimate_model(formula=formula, data=model_data, A=model_data$X, spde=spde, prec_initial=1)

      mu.tmp <- INLA_result$misc$theta.mode #for joint group model
      Q.tmp <- solve(INLA_result$misc$cov.intern) #for joint group model

      mu.theta.tmp <- mu.theta.tmp + as.vector(Q.tmp%*%mu.tmp)
      Q.theta <- Q.theta + Q.tmp
      theta.sub <- cbind(theta.sub, res.hyper$mode)
      rm(mu.tmp, Q.tmp)
    }
    mu.theta <- solve(Q.theta, mu.theta.tmp)

  # Drawing samples from q(theta|y)

  nsamp <- 50
  logwt <- rep(NA, nsamp)

  theta.tmp <- mvrnorm(nsamp, mu.theta, solve(Q.theta))
  for(i in 1:nsamp){ logwt[i] <- F.logwt(theta.tmp[i,], spde[[h]], mu.theta, Q.theta, M) }

  #weights to apply to each posterior sample of theta
  wt.tmp <- exp(logwt - max(logwt))
  wt <- wt.tmp/(sum(wt.tmp))


  ## Create index vectors
  n.mesh <- mesh$n
  ind_beta <- list()
  for(k in 1:K){
    ind_beta[[k]] <- 1:n.mesh + (k-1)*n.mesh
  }

  #get posterior quantities of beta, conditional on a value of theta
  if(is.null(no_cores)) {
    parallel <- FALSE
    beta.post.samps <- apply(theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, K=K, M, Xcros.all, Xycros.all, thresholds=thresholds, alpha=0.01, ind_beta=ind_beta)
  } else {
    max_no_cores <- min(detectCores() - 1, 25)
    no_cores <- min(max_no_cores, no_cores)
    cl <- makeCluster(no_cores)
    beta.post.samps <- parApply(cl, theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, K=K, M, Xcros.all, Xycros.all, thresholds=thresholds, alpha=0.01, ind_beta=ind_beta)
    print(Sys.time() - t0)
    stopCluster(cl)
  }
  #in sequence, 8 min per iteration for motor task,  5-6 min for gambling task (for 20 subjects and 3 activation thresholds)
  #in parallel, 21 min total for gambling task!
  #in parallel, 24 min for motor task!
  #with 50 iterations, we save 50*8 - 25 = 375 min (6.25 hours!)


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
  probs.all <- array(0, dim=c(n.mesh, K, U)) #last dimension is for different activation thresholds

  #posterior mean
  beta.pop <- as.vector(mu.tot%*%wt)

  for(k in 1:K){
      beta.pop.k <- beta.pop[ind_beta[[k]]]
      betas.all[,k] <- as.vector(beta.pop.k)
  }

  #posterior probabilities
  for(u in 1:U){
    for(k in 1:K){
      F.pop.uk <- as.vector(F.tot[[u]][[k]]%*%wt)
      probs.all[,k,u] <- as.vector(F.pop.uk)
    }
  }

  result <- list(beta_estimates = betas.all, active_beta_probs = probs.all)

  class(result) <- "BayesGLM_group"

  return(result)

}#end loop over hemispheres


