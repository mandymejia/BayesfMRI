#'  BayesGLMEM for 2D slice
#'
#'  Spatial Bayesian GLM expectation maximization for fMRI task activation on
#'    2d slice volumetric data.
#'
#' @param BOLD A list of sessions, each with a three-dimensional array in which
#'   the first two dimensions correspond to the size of the fMRI slice in space
#'   and the last dimension corresponds to time
#' @param binary_mask (optional) a binary brain slice image used to mask
#'   the BOLD data and make a more efficient network mesh for the
#'   neighborhood definitions
#' @param design,onsets,TR All or none must be provided.
#'
#'   \code{design} is a \eqn{T x K} task design matrix (or list of such
#'   matrices, for multiple-session modeling) with column names representing
#'   tasks. Each column represents the expected BOLD response due to each task,
#'   a convolution of the hemodynamic response function (HRF) and the task
#'   stimulus. Note that the scale of the regressors will affect the scale and
#'   interpretation of the beta coefficients, so imposing a proper scale (e.g.,
#'   set maximum to 1) is recommended.
#'
#'   \code{onsets} is a matrix of onsets (first column) and durations (second column)
#'   for each task in seconds, organized as a list where each element of the
#'   list corresponds to one task. Names of list should be task names. (Or for
#'   multi-session modeling, a list of such lists.)
#'
#'   \code{TR} is the temporal resolution of the data in seconds.
#' @param nuisance (Optional) A TxJ matrix of nuisance signals (or list of such
#'   matrices, for multiple-session modeling).
#' @param nuisance_include (Optional) Additional nuisance covariates to include.
#'   Default is 'drift' (linear and quadratic drift terms) and 'dHRF' (temporal
#'   derivative of each column of design matrix).
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @inheritParams num.threads_Param
#' @param GLM_method Either 'Bayesian' for spatial Bayesian GLM only,
#'   'classical' for the classical GLM only, or 'both' to return both classical
#'   and Bayesian estimates of task activation.
#' @param session_names (Optional) A vector of names corresponding to each
#'   session.
#' @param outfile (Optional) File name (without extension) of output file for
#'   BayesGLMEM result to use in Bayesian group modeling.
#' @inheritParams verbose_Param_direct_TRUE
#'
#' @importFrom utils head
#' @importFrom INLA inla.spde.make.A
#'
#' @return An object of class \code{"BayesGLM"}, a list containing...
#'
#' @export
BayesGLMEM_slice <- function(
  BOLD,
  design = NULL,
  onsets=NULL,
  TR=NULL,
  nuisance=NULL,
  nuisance_include=c('drift','dHRF'),
  binary_mask = NULL,
  scale_BOLD = TRUE,
  scale_design = TRUE,
  num.threads = 4,
  tol = 1e-4,
  GLM_method = 'both',
  session_names = NULL,
  outfile = NULL,
  verbose = FALSE) {

  do_Bayesian <- (GLM_method %in% c('both','Bayesian'))
  do_classical <- (GLM_method %in% c('both','classical'))

  image_dims <- head(dim(BOLD[[1]]),-1)
  if (is.null(binary_mask))
    binary_mask <- matrix(1, nrow = image_dims[1], ncol = image_dims[2])

  mesh <- make_slice_mesh(binary_mask)
  # Create a conversion matrix
  in_binary_mask <- which(binary_mask == 1, arr.ind = T)
  in_binary_mask <- in_binary_mask[,2:1]
  convert_mat_A <- INLA::inla.spde.make.A(mesh = mesh, loc = in_binary_mask)

  # Name sessions and check compatibility of multi-session arguments
  n_sess <- length(BOLD)
  if(n_sess==1){
    if(is.null(session_names)) session_names <- 'single_session'
  } else {
    if(is.null(session_names)) session_names <- paste0('session', 1:n_sess)
  }
  if(length(session_names) != n_sess)
    stop('If session_names is provided, it must be of the same length as BOLD')

  cat('\n SETTING UP DATA \n')

  if(is.null(design)) {
    make_design <- TRUE
    design <- vector('list', length=n_sess)
  } else {
    make_design <- FALSE
  }

  for(ss in 1:n_sess){
    if(make_design){
      cat(paste0('    Constructing design matrix for session ', ss, '\n'))
      design[[ss]] <- make_HRFs(onsets[[ss]], TR=TR, duration=ntime)
    }
  }

  ### Check that design matrix names consistent across sessions
  if(n_sess > 1){
    tmp <- sapply(design, colnames)
    tmp <- apply(tmp, 1, function(x) length(unique(x)))
    if(max(tmp) > 1)
      stop('task names must match across sessions for multi-session modeling')
  }

  cat('\n RUNNING MODEL \n')

  classicalGLM <- NULL
  BayesGLM <- NULL

  ### FORMAT DESIGN MATRIX
  for(ss in 1:n_sess){
    if(scale_design){
      design[[ss]] <- scale_design_mat(design[[ss]])
    } else {
      design[[ss]] <- scale(design[[ss]], scale=FALSE) #center design matrix
      # to eliminate baseline
    }
  }

  ### ADD ADDITIONAL NUISANCE REGRESSORS
  if(!is.null(nuisance)) {
    for (ss in 1:n_sess) {
      ntime <- nrow(design[[ss]])
      if ('drift' %in% nuisance_include) {
        drift <- (1:ntime) / ntime
        nuisance[[ss]] <-
          cbind(nuisance[[ss]], drift, drift ^ 2)
      }
      if ('dHRF' %in% nuisance_include) {
        dHRF <- gradient(design[[ss]])
        nuisance[[ss]] <-
          cbind(nuisance[[ss]], dHRF)
      }
    }
  } else {
    nuisance <- list()
    for (ss in 1:n_sess) {
      ntime <- nrow(design[[ss]])
      if ('drift' %in% nuisance_include) {
        drift <- (1:ntime) / ntime
        nuisance[[ss]] <- cbind(drift, drift ^ 2)
      }
      if ('dHRF' %in% nuisance_include) {
        dHRF <- gradient(design[[ss]])
        nuisance[[ss]] <- dHRF
      }
    }
  }

  scale_design <- F # This is done to prevent double-scaling in BayesGLMEM

  #set up session list
  session_data <- vector('list', n_sess)
  names(session_data) <- session_names
  for(ss in 1:n_sess){
    # BOLD[[ss]] <- scale_timeseries(BOLD = BOLD[[ss]],scale = scale_BOLD,transpose = TRUE)
    # sess <- list(BOLD = as.matrix(BOLD[[ss]]%*%convert_mat_A), design=design[[ss]])
    sess <- list(BOLD = BOLD[[ss]], design=design[[ss]])
    if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
    session_data[[ss]] <- sess
  }

  cat(str(binary_mask),"\n")

  ### FIT GLM(s)

  if(do_classical) classicalGLM_out <- classicalGLM(session_data,
                                                    scale_BOLD=scale_BOLD,
                                                    scale_design = scale_design)
  if(do_Bayesian) BayesGLM_out <- BayesGLMEM(session_data,
                                           mesh = mesh,
                                           mask = binary_mask,
                                           scale_BOLD=scale_BOLD,
                                           scale_design = scale_design,
                                           num.threads=num.threads,
                                           tol = tol,
                                           outfile = outfile,
                                           verbose=verbose)


  # Extract the point estimates
  cat(str(convert_mat_A),"\n")
  cat(str(BayesGLM_out$beta_estimates),"\n")
  point_estimates <- sapply(session_names, function(sn){
    as.matrix(convert_mat_A %*% BayesGLM_out$beta_estimates[[sn]])
  }, simplify = F)

  classical_slice <- Bayes_slice <- vector('list', n_sess)
  names(classical_slice) <- names(Bayes_slice) <- session_names
  for(ss in 1:n_sess){
    num_tasks <- ncol(design[[ss]])
    if(do_classical){
      classical_slice[[ss]] <- sapply(seq(num_tasks), function(tn) {
        image_coef <- binary_mask
        not_na <- which(!is.na(classicalGLM_out[[ss]][,tn]))
        image_coef[image_coef == 1] <- classicalGLM_out[[ss]][not_na,tn]
        image_coef[binary_mask == 0] <- NA
        return(image_coef)
      },simplify = F)
    }
    if(do_Bayesian){
      # mat_coefs <- as.matrix(convert_mat_A %*% BayesGLM_out$beta_estimates[[ss]])
      Bayes_slice[[ss]] <- sapply(seq(num_tasks), function(tn) {
        image_coef <- binary_mask
        image_coef[image_coef == 1] <- point_estimates[[ss]][,tn]
        image_coef[binary_mask == 0] <- NA
        return(image_coef)
      },simplify = F)
      names(Bayes_slice[[ss]]) <- beta_names
    }
  }

  if (do_Bayesian) {
    beta_names <- BayesGLM_out$beta_names
  } else {
    beta_names <- NULL
  }

  result <- list(session_names = session_names,
                 beta_names = beta_names,
                 betas_Bayesian = Bayes_slice,
                 betas_classical = classical_slice,
                 GLMs_Bayesian = BayesGLM_out,
                 GLMs_classical = classicalGLM_out,
                 design = design,
                 mask = binary_mask)
  class(result) <- "BayesGLM_slice"
  return(result)
}

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
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#' @inheritParams mesh_Param_inla
#' @param mask (Optional) A length \eqn{V} logical vector indicating if each
#'  vertex is to be included.
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @param num.threads Not used yet.
#' @param outfile File name where results will be written (for use by
#'  \code{BayesGLM2}).
#' @inheritParams verbose_Param_direct_TRUE
#'
#' @return A list containing...
#'
#' @importFrom INLA inla.spde2.matern
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#'
#' @export
BayesGLMEM <- function(data,
                       vertices = NULL,
                       faces = NULL,
                       mesh = NULL,
                       mask = NULL,
                       scale_BOLD = TRUE,
                       scale_design = TRUE,
                       num.threads = 4,
                       tol = 1e-4,
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
  K <- ncol(data[[1]]$design) #number of tasks
  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
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
  in_mask <- which(mask == 1, arr.ind = T)
  in_mask <- in_mask[,2:1]
  spde <- INLA::inla.spde2.matern(mesh)

  #collect data and design matrices
  y_all <- c()
  X_all_list <- NULL
  design <- vector('list', length=n_sess)

  for(s in 1:n_sess){

    #extract and mask BOLD data for current session
    BOLD_s <- data[[s]]$BOLD

    #scale data to represent % signal change (or just center if scale=FALSE)
    BOLD_s <- scale_timeseries(BOLD_s, scale=scale_BOLD, transpose = TRUE)
    if(scale_design) {
      design_s <- scale_design_mat(data[[s]]$design)
    } else {
      design_s <- scale(data[[s]]$design, scale = F)
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
    data_org <- organize_data(y_reg, X_reg)
    y_vec <- data_org$y
    X_list <- list(data_org$X)
    names(X_list) <- session_names[s]

    y_all <- c(y_all, y_vec)
    X_all_list <- c(X_all_list, X_list)
  }

  #construct betas and repls objects
  replicates_list <- BayesfMRI:::organize_replicates(n_sess=n_sess, n_task=K, mesh=mesh)
  betas <- replicates_list$betas
  repls <- replicates_list$repls

  #organize the formula and data objects
  beta_names <- names(betas)
  repl_names <- names(repls)
  n_beta <- length(names(betas))
  model_data <- BayesfMRI:::make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)
  # > Model setup and initials ----
  Psi_k <- INLA::inla.spde.make.A(mesh = mesh, loc = in_mask)
  Psi <- Matrix::bdiag(rep(list(Psi_k),K))
  A <- Matrix::crossprod(model_data$X%*%Psi)
  # Initial values for kappa and tau
  kappa2_new <- rep(4,K) # This is a value that matches BayesGLM
  phi_new <- 1 / (4*pi*kappa2_new*4) # This is a value that matches BayesGLM
  sigma2_new <- var(model_data$y)
  Q_k <- mapply(spde_Q_phi, kappa2 = kappa2_new, phi = phi_new, MoreArgs = list(spde = spde), SIMPLIFY = F)
  Q_new <- Matrix::bdiag(Q_k)
  Sig_inv <- Q_new + A/sigma2_new
  m <- Matrix::crossprod(model_data$X%*%Psi,model_data$y) / sigma2_new
  mu <- Matrix::solve(Sig_inv,m)
  Q_chol <- Matrix::chol(Q_new, pivot = T)
  Q_chol_inv <- Matrix::solve(Q_chol)
  Sigma_new <- Q_chol_inv %*% (A/sigma2_new + diag(1,nrow(A))) %*% t(Q_chol_inv)
  # mean model log-likelihood (for checking)
  llik_new <- -length(model_data$y)/2 * log(sigma2_new) -
    crossprod(model_data$y - mean(model_data$y))/(2*sigma2_new)
  if(verbose) {
    cat("Log-likelihood for the mean model:", llik_new,"\n")
  }
  llik <- -Inf
  # > EM algorithm ----
  cat(".... performing the EM algorithm\n")
  # diff_Theta <- Inf
  step <- 1
  n <- spde$n.spde
  while(llik_new > llik + tol) {
    Q <- Q_new
    kappa2 <- kappa2_new
    phi <- phi_new
    sigma2 <- sigma2_new
    Sigma <- Sigma_new
    llik <- llik_new
    for(k in seq(K)) {
      # k_inds <- seq(V) + (k-1)*V
      k_inds <- seq(n) + (k-1)*n
      Qp <- Q_prime(kappa2[k], spde)
      Tr_QEww <- (sum(Matrix::colSums(Qp%*%Sigma[k_inds,k_inds])) +
                    crossprod(mu[k_inds,],Q[k_inds,k_inds])%*%mu[k_inds,])@x
      phi_new[k] <- Tr_QEww / (4*pi*V)
      optim_output_k <-
        optim(
          par = kappa2[k],
          fn = neg_kappa_fn,
          method = "L-BFGS-B",
          spde = spde,
          phi = phi[k],
          Sigma = Sigma[k_inds,k_inds],
          mu = mu[k_inds,],
          lower = 1e-4
        )
      kappa2_new[k] <- optim_output_k$par
    }
    Qk_new <- mapply(spde_Q_phi,kappa2 = kappa2_new, phi = phi_new,
                     MoreArgs = list(spde=spde), SIMPLIFY = F)
    Q_new <- Matrix::bdiag(Qk_new)

    Sig_inv <- Q_new + A/sigma2
    m <- t(model_data$X%*%Psi)%*%model_data$y / sigma2
    mu <- solve(Sig_inv,m)
    Q_chol <- chol(Q, pivot = T)
    Q_chol_inv <- solve(Q_chol)
    Sigma <- Q_chol_inv %*% (A/sigma2 + diag(1,nrow(A))) %*% t(Q_chol_inv)


    sigma2_new <-
      as.numeric(crossprod(model_data$y) -
                   2*crossprod(model_data$y,model_data$X%*%Psi%*%mu) +
                   crossprod(model_data$X%*%Psi%*%mu) +
                   sum(colSums(A*Sigma))) / length(model_data$y)

    llik_new <- -length(model_data$y)/2 *  log(sigma2) -
      crossprod(model_data$y - as.numeric(model_data$X%*%Psi%*%mu))/(2*sigma2)
    if(verbose) {
      cat("Step",step, "kappa^2 =",kappa2_new, "phi =", phi_new, "sigma^2 =",sigma2_new ,
          "log-likelihood =",llik_new,"\n")
    }
    step <- step+1
  }
  cat(".... EM algorithm complete!")

  Sig_inv <- Q + A/sigma2
  m <- t(model_data$X%*%Psi)%*%model_data$y / sigma2
  mu <- solve(Sig_inv,m)
  Q_chol <- chol(Q, pivot = T)
  Q_chol_inv <- solve(Q_chol)
  Sigma_new <- Q_chol_inv %*% (A/sigma2 + diag(1,nrow(A))) %*% t(Q_chol_inv)

  beta_estimates <- matrix(mu@x,nrow = n, ncol = K)
  colnames(beta_estimates) <- beta_names
  beta_estimates <- list(beta_estimates)
  names(beta_estimates) <- session_names
  theta_estimates <- c(sigma2,c(phi,kappa2))
  names(theta_estimates) <- c("sigma2",paste0("phi_",seq(K)),paste0("kappa2_",seq(K)))
  #extract stuff needed for group analysis
  tau2 <- 1 / (4*pi*kappa2*phi)
  mu.theta <- c(c(rbind(tau2,kappa2)), sigma2) # This is a guess about the order and might be wrong
  Q.theta <- Q

  # Construct object to be returned
  result <- list(mesh = mesh,
                 mask = mask,
                 design = design,
                 session_names = session_names,
                 beta_names = beta_names,
                 beta_estimates = beta_estimates,
                 theta_estimates = theta_estimates,
                 mu.theta = mu.theta, #for joint group model
                 Q.theta = Q.theta, #for joint group model
                 y = y_all, #for joint group model
                 X = X_all_list, #for joint group model
                 call = match.call())

  class(result) <- "BayesGLM"

  if(!is.null(outfile)){
    saveRDS(result, file=outfile)
  }

  return(result)
}
