#' Applies joint approach to group-level analysis to task fMRI data 
#'
#' @param data A list of signle-session data of M subjects, where each subject is a list with elements
#' BOLD, design and nuisance. See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @param vertices A Vx3 matrix of vertex locations of the triangular mesh in Euclidean space.
#' @param faces A Wx3 matrix, where each row contains the vertex indices for a given face or triangle in the triangular mesh.
#' @param mesh A `inla.mesh` object.  Must be provided if and only if `vertices` and `faces` are not.
#' @param mask Sarah knows what this does!
#' @param scale If TRUE, scale timeseries data so estimates represent percent signal change.  Else, do not scale.
#' @param thresholds A vector of thresholds for activation maps.
#' @param method The type of methods: joint and ...
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern
#'
#' @examples \dontrun{}
BayesGLM_group <- function(data, method = NULL, vertices = NULL, faces = NULL, mesh = NULL, mask = NULL, scale=TRUE,  thresholds=c(0,0.5,1)){
  #check whether data is a list OR a session (for single-session analysis)
  #check whether each element of data is a session (use is.session)
  # V = number of data locations
  # T = length of time series for each session (vector)
  # K = number of unique tasks in all sessions

  #need to check that sessions are consistent in terms of V, K?

  #INLA:::inla.dynload.workaround() #avoid error on creating mesh

  # Check to see that the INLA package is installed
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("This function requires the INLA package (see www.r-inla.org/download)")


  # Check to see if PARDISO is installed
  if(!exists("inla.pardiso.check", mode = "function")){
    warning("Please update to the latest version of INLA for full functionality and PARDISO compatibility (see www.r-inla.org/download)")
  }else{
  if(inla.pardiso.check() == "FAILURE: PARDISO IS NOT INSTALLED OR NOT WORKING"){
  warning("Consider enabling PARDISO for faster computation (see inla.pardiso())")}
  #inla.pardiso()
 }


  #check that only mesh OR vertices+faces supplied
  has_mesh <- !is.null(mesh)
  has_verts_faces <- !is.null(vertices) & !is.null(faces)
  has_howmany <- has_mesh + has_verts_faces
  if(has_howmany != 1) stop('Must supply EITHER mesh OR vertices and faces.')

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks

  if(!is.list(data)) stop('I expect data to be a list, but it is not')
    data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

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
  no_cores <- min(detectCores() - 1, 25)
  cl <- makeCluster(no_cores)
  t0 <- Sys.time()
  #in sequence, 8 min per iteration for motor task,  5-6 min for gambling task (for 20 subjects and 3 activation thresholds)
  #in parallel, 21 min total for gambling task!
  #in parallel, 24 min for motor task!
  #with 50 iterations, we save 50*8 - 25 = 375 min (6.25 hours!)
  U <- length(thresholds)
  beta.post.samps <- parApply(cl, theta.tmp, MARGIN=1, FUN=beta.posterior.thetasamp, spde=spde, K=K, M, Xcros.all, Xycros.all, thresholds=thresholds, alpha=0.01, ind_beta=ind_beta)
  print(Sys.time() - t0)
  stopCluster(cl)


  #organize samples
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


