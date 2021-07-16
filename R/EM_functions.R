#' Expectation Maximization for the Bayesian General Linear Model
#'
#' Applies spatial Bayesian GLM to task fMRI data
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param data A list of sessions, where each session is a list with elements
#'  BOLD, design and nuisance. See \code{?create.session} and \code{?is.session}
#'  for more details.
#' List element names represent session names.
#' @param beta_names (Optional) Names of tasks represented in design matrix
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#' @inheritParams mesh_Param_inla
#' @param mask (Optional) A length \eqn{V} logical vector indicating if each
#'  vertex is to be included.
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @param EM_method Either "joint" or "separate" for choosing whether covariates
#'   should share hyperparameter values.
#' @param use_SQUAREM (logical) Should the SQUAREM package be used to speed up
#'   convergence?
#' @param tol If use_SQUAREM == TRUE, an absolute change limit for
#'   when the EM algorithm should be stopped (Default = 1e-3). If use_SQUAREM ==
#'   FALSE, a percent change limit for when the EM algorithm should be stopped
#'   (Default = 1). A value of NULL will result in the default value being used.
#' @param num.threads (optional) allows users to specify the number of threads used
#'   to work in parallel across the different tasks
#' @param outfile (Optional) File name (without extension) of output file for
#'   BayesGLMEM result to use in Bayesian group modeling.
#' @inheritParams verbose_Param_direct_TRUE
#'
#' @return A list containing...
#'
#' @importFrom INLA inla.spde2.matern inla.qinv inla.qsolve
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#' @importFrom SQUAREM squarem
#' @importFrom parallel makeCluster parApply
#'
#' @export
BayesGLMEM <- function(data,
                       beta_names = NULL,
                       vertices = NULL,
                       faces = NULL,
                       mesh = NULL,
                       mask = NULL,
                       scale_BOLD = TRUE,
                       scale_design = TRUE,
                       EM_method = "separate",
                       use_SQUAREM = TRUE,
                       tol = NULL,
                       num.threads = 1,
                       outfile = NULL,
                       verbose = FALSE) {

  # > Data setup ----
  #check whether data is a list OR a session (for single-session analysis)
  #check whether each element of data is a session (use is.session)
  # V = number of data locations
  # T = length of time series for each session (vector)
  # K = number of unique tasks in all sessions

  #check that only mesh OR vertices+faces supplied
  has_mesh <- !is.null(mesh)
  has_verts_faces <- !is.null(vertices) & !is.null(faces)
  has_howmany <- has_mesh + has_verts_faces
  if(has_howmany != 1) stop('Must supply EITHER mesh OR vertices and faces.')

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data)
  n_sess <- length(session_names)

  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  V <- ncol(data[[1]]$BOLD) #number of data locations
  is_missing <- is.na(data[[1]]$BOLD[1,])
  V_nm <- V - sum(is_missing)
  ntime <- nrow(data[[1]]$BOLD)
  is_pw <- nrow(data[[1]]$design) == (ntime * sum(!is_missing))
  if(is_pw) {
    K <- ncol(data[[1]]$design) / V_nm # number of tasks
  } else {
    K <- ncol(data[[1]]$design) #number of tasks
  }
  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(is_pw) {
      if(ncol(data[[s]]$design) / V_nm != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
    } else {
      if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
    }

  }

  if(!is.null(beta_names)){
    if(length(beta_names) != K) stop(paste0('I detect ', K, ' task based on the design matrix, but the length of beta_names is ', length(beta_names), '.  Please fix beta_names.'))
  }

  if(is.null(beta_names)){
    if(!is_pw){
      beta_names_maybe <- colnames(data[[1]]$design) #if no prewhitening, can grab beta names from design (if not provided)
      if(!is.null(beta_names_maybe)) beta_names <- beta_names_maybe
      if(is.null(beta_names_maybe)) beta_names <- paste0('beta',1:K)
    } else {
      beta_names <- paste0('beta',1:K)
    }
  }

  if(is.null(outfile)){
    message('No value supplied for `outfile`, which is required for post-hoc group modeling.')
  }

  if(is.null(mesh)) mesh <- make_mesh(vertices, faces)

  #ID any zero-variance voxels and remove from analysis
  zero_var <- sapply(data, function(x){
    x$BOLD[is.na(x$BOLD)] <- 0 #to detect medial wall locations coded as NA
    x$BOLD[is.nan(x$BOLD)] <- 0 #to detect medial wall locations coded as NaN
    vars <- matrixStats::colVars(x$BOLD)
    return(vars < 1e-6)
  })
  zero_var <- (rowSums(zero_var) > 0) #check whether any vertices have zero variance in any session

  #1. Apply mask to mesh, data and zero_var
  #2. If sum(zero_var) > 0, remove zero_var locations from data and create Amat
  #   Else, let Amat = identity matrix

  if(is.null(mesh) | is.null(mask)) {
    if(sum(zero_var) > 0){
      if(!is.null(mask)) mask[zero_var==TRUE] <- 0
      if(is.null(mask)) mask <- !zero_var
    }

    if(!is.null(mask) & sum(mask) != V) {
      mask <- as.logical(mask)
      mesh <- excursions::submesh.mesh(mask, mesh)
      mesh$idx$loc <- mesh$idx$loc[!is.na(mesh$idx$loc)]
      for(s in 1:n_sess){
        data[[s]]$BOLD <- data[[s]]$BOLD[,mask]
      }
      V <- sum(mask)
    }
  }
  in_mask <- NULL
  if("matrix" %in% class(mask)) {
    in_mask <- which(mask == 1, arr.ind = T)
    in_mask <- in_mask[,2:1]
  }
  spde <- INLA::inla.spde2.matern(mesh)

  #collect data and design matrices
  y_all <- c()
  X_all_list <- NULL
  design <- vector('list', length=n_sess)

  for(s in 1:n_sess){

    #extract and mask BOLD data for current session
    BOLD_s <- data[[s]]$BOLD

    #scale data to represent % signal change (or just center if scale=FALSE)
    if(!is_pw) BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale_BOLD)
    if(scale_design) {
      design_s <- scale_design_mat(data[[s]]$design)
    } else {
      if(!is_pw) {
        design_s <- scale(data[[s]]$design, scale = F)
      } else {
        design_s <- data[[s]]$design # Don't scale prewhitened data or the matrix will not be sparse
      }
    }
    design[[s]] <- design_s #after scaling but before nuisance regression

    #regress nuisance parameters from BOLD data and design matrix
    if('nuisance' %in% names(data[[s]])){
      nuisance_s <- data[[s]]$nuisance
      y_reg <- nuisance_regress(BOLD_s, nuisance_s)
      X_reg <- nuisance_regress(design_s, nuisance_s)
    } else {
      y_reg <- BOLD_s
      X_reg <- design_s
    }

    #set up data and design matrix
    if(!is_pw) {
      data_org <- organize_data(y_reg, X_reg)
    } else {
      data_org <- organize_data_pw(y_reg, X_reg)
    }
    y_vec <- data_org$y
    X_list <- list(data_org$X)
    names(X_list) <- session_names[s]

    y_all <- c(y_all, y_vec)
    X_all_list <- c(X_all_list, X_list)
  }

  #construct betas and repls objects
  replicates_list <- organize_replicates(n_sess=n_sess, beta_names = beta_names, mesh=mesh)
  betas <- replicates_list$betas
  repls <- replicates_list$repls

  #organize the formula and data objects
  beta_names <- names(betas)
  repl_names <- names(repls)
  n_beta <- length(names(betas))
  model_data <- make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)
  # > Model setup and initials ----
  if(is.null(tol)){
    if(use_SQUAREM) tol <- 1e-3
    if(!use_SQUAREM) tol <- 1
  }
  if(!is.null(in_mask)) {
    Psi_k <- INLA::inla.spde.make.A(mesh = mesh, loc = in_mask)
    Psi <- Matrix::bdiag(rep(list(Psi_k),K))
    A <- Matrix::crossprod(model_data$X%*%Psi)
  } else {
    Psi <- Matrix::Diagonal(n = ncol(model_data$X))
    A <- Matrix::crossprod(model_data$X)
  }
  # Initial values for kappa and tau
  # Using values matching BayesGLM
  if(EM_method == "joint") num.threads <- 1
  kappa2 <- 4
  phi <- 1 / (4*pi*kappa2*4) # This is a value that matches BayesGLM
  # sigma2 <- var(model_data$y)
  # theta <- c(kappa2, phi, sigma2)
  # Using values based on the classical GLM
  cat("... FINDING BEST GUESS INITIAL VALUES\n")
  beta_hat <- (Matrix::solve(Matrix::crossprod(model_data$X)) %*%
                 Matrix::crossprod(model_data$X,model_data$y))@x
  sigma2 <- ((Matrix::crossprod(model_data$y) -
                Matrix::crossprod(model_data$y,model_data$X) %*% Psi %*%
                beta_hat + sum(diag(Matrix::crossprod(model_data$X%*%Psi) %*%
                                      Matrix::tcrossprod(beta_hat)))) /
               length(model_data$y))@x
  if(EM_method == "joint") {
    # require(SQUAREM)
    init_output <-
      SQUAREM::squarem(
        par = c(kappa2, phi),
        fixptfn = init_fixpt,
        # objfn = init_objfn, # This isn't strictly necessary, and may cost a small amount of time.
        spde = spde,
        beta_hat = beta_hat,
        control = list(tol = 1e-3, trace = verbose)
      )
    theta <- c(init_output$par, sigma2)
    cat("...... DONE!\n")
  }
  if(EM_method == "separate") {
    beta_hat <- matrix(beta_hat, ncol = K)
    # require(parallel)
    cl <- parallel::makeCluster(min(num.threads,K))
    kappa2_phi <- parallel::parApply(cl,beta_hat,2, function(bh, kappa2, phi, spde, verbose) {
      # require(SQUAREM)
      # require(Matrix)
      init_output <-
        SQUAREM::squarem(
          par = c(kappa2, phi),
          fixptfn = init_fixpt,
          # objfn = init_objfn, # Not needed
          spde = spde,
          beta_hat = bh,
          control = list(tol = 1e-3, trace = verbose)
        )
      return(init_output)
    },kappa2 = kappa2, phi = phi, spde = spde, verbose = verbose)
    kappa2_phi <- sapply(kappa2_phi,function(x) x$par)
    theta <- c(t(kappa2_phi),sigma2)
    cat("...... DONE!\n")
    parallel::stopCluster(cl)
  }

  # > Start EM algorithm ----
  if(EM_method == "joint") {
    if(length(theta) != 3) stop("The length of theta should be 3 for the joint update")
    em_fn <- GLMEM_fixptjoint
    k_idx <- 1
    p_idx <- 2
    s_idx <- 3
  }
  if(EM_method == "separate") {
    em_fn <- GLMEM_fixptseparate
    k_idx <- seq(K)
    p_idx <- seq(K) + K
    s_idx <- (2*K + 1)
  }
  if(use_SQUAREM) {
    squareem_output <-
      squarem(
        par = theta,
        fixptfn = em_fn,
        # objfn = GLMEM_objfn,
        control = list(tol = tol, trace = verbose),
        spde = spde,
        model_data = model_data,
        Psi = Psi,
        K = K,
        A = A,
        num.threads = num.threads
      )
    theta_new <- squareem_output$par
    em_output <- list(
      theta_new = theta_new,
      kappa2_new = theta_new[k_idx],
      phi_new = theta_new[p_idx],
      sigma2_new = theta_new[s_idx]
    )
  } else {
    step <- 1
    max_pct_change <- Inf
    while(max_pct_change > tol | step <= 5) {
      theta_new <-
        em_fn(
          theta = theta,
          spde = spde,
          model_data = model_data,
          Psi = Psi,
          K = K,
          A = A,
          num.threads = num.threads
        )
      kappa2_new <- theta_new[k_idx]
      phi_new <- theta_new[p_idx]
      sigma2_new <- theta_new[s_idx]
      sigma2_pct_change <- 100*abs((sigma2_new - sigma2) / sigma2)
      phi_pct_change <- 100*abs((phi_new - phi) / phi)
      kappa2_pct_change <- 100*abs((kappa2_new - kappa2) / kappa2)
      max_pct_change <- max(sigma2_pct_change,phi_pct_change,kappa2_pct_change)
      if(verbose) {
        cat(
          "Step",
          step,
          "kappa^2 (%change) =",
          kappa2_new,
          "(",
          kappa2_pct_change,
          ") phi (%change) =",
          phi_new,
          "(",
          phi_pct_change,
          ") sigma^2 (%change) =",
          sigma2_new,
          "(",
          sigma2_pct_change,
          ")",
          "\n"
        )
      }
      kappa2 <- kappa2_new
      phi <- phi_new
      sigma2 <- sigma2_new
      theta <- theta_new
      step <- step+1
    }
    em_output <- list(
      theta_new = theta_new,
      kappa2_new = kappa2_new,
      phi_new = phi_new,
      sigma2_new = sigma2_new
    )
  }
  # > End EM algorithm ----
  cat(".... EM algorithm complete!\n")
  list2env(em_output, envir = environment())
  Qk_new <- mapply(spde_Q_phi,kappa2 = kappa2_new, phi = phi_new,
                   MoreArgs = list(spde=spde), SIMPLIFY = F)
  if(EM_method == "joint") {
    Q <- Matrix::bdiag(rep(Qk_new,K))
  } else {
    Q <- Matrix::bdiag(Qk_new)
  }
  if(n_sess > 1) Q <- Matrix::bdiag(lapply(seq(n_sess), function(x) Q))
  Sig_inv <- Q + A/sigma2_new
  m <- Matrix::t(model_data$X%*%Psi)%*%model_data$y / sigma2_new
  mu <- INLA::inla.qsolve(Sig_inv,m)
  Sigma_new <- INLA::inla.qsolve(Sig_inv, Matrix::Diagonal(n = nrow(Sig_inv)), method = "solve")

  beta_estimates <- matrix(mu,nrow = spde$n.spde, ncol = K*n_sess)
  colnames(beta_estimates) <- rep(beta_names, n_sess)
  beta_estimates <- lapply(seq(n_sess), function(ns) beta_estimates[,(seq(K) + K * (ns - 1))])
  names(beta_estimates) <- session_names
  theta_estimates <- c(sigma2_new,c(phi_new,kappa2_new))
  if(EM_method == "joint") {
    names(theta_estimates) <- c("sigma2","phi","kappa2")
  } else {
    names(theta_estimates) <- c("sigma2",paste0("phi_",seq(K)),paste0("kappa2_",seq(K)))
  }
  #extract stuff needed for group analysis
  tau2 <- 1 / (4*pi*kappa2_new*phi_new)
  mu.theta <- c(log(sigma2_new),c(rbind(log(sqrt(tau2)),log(sqrt(kappa2_new))))) # This is a guess about the order and might be wrong
  # Q.theta <- Q # This is not right. This is supposed to be the covariance between
  # the hyperparameters (kappa,phi,sigma2) This might need to be examined. Perhaps
  # an estimate can be made using the iteration values for the parameters?

  # Construct object to be returned
  result <- list(mesh = mesh,
                 mask = mask,
                 design = design,
                 session_names = session_names,
                 beta_names = beta_names,
                 beta_estimates = beta_estimates,
                 theta_estimates = theta_estimates,
                 posterior_Sig_inv = Sig_inv, # For excursions
                 mu.theta = mu.theta, #for joint group model
                 y = y_all, #for joint group model
                 X = X_all_list, #for joint group model
                 call = match.call())
  class(result) <- "BayesGLM"
  return(result)
}

