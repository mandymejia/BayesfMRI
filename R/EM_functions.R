#' Applies spatial Bayesian GLM to task fMRI data using the EM algorithm
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @param vertices A Vx3 matrix of vertex locations of the triangular mesh in Euclidean space.
#' @param faces A Wx3 matrix, where each row contains the vertex indices for a given face or triangle in the triangular mesh. W is the number of faces in the mesh.
#' @param mesh A `inla.mesh` object.  Must be provided if and only if `vertices` and `faces` are not.
#' @param mask (Optional) A logical or 0/1 vector of length V indicating which vertices are to be included.
#' @param scale_BOLD If TRUE, scale timeseries data so estimates represent percent signal change.  Else, center but do not scale.
#' @param scale_design If TRUE, scale the design matrix by dividing each column
#'   by its maximum value, and then subtracting the new column mean.
#' @param num.threads Maximum number of threads the inla-program will use for model estimation
#' @param return_INLA_result If TRUE, object returned will include the INLA model object (can be large).  Default is FALSE.  Required for running \code{id_activations} on \code{BayesGLM} model object (but not for running \code{BayesGLM_joint} to get posterior quantities of group means or contrasts).
#' @param outfile File name where results will be written (for use by \code{BayesGLM_group}).
#' @param verbose Logical indicating if INLA should run in a verbose mode (default FALSE).
#' @param contrasts A list of contrast vectors to be passed to
#'   \code{\link{inla}}.
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern inla.pardiso.check inla.setOption
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
BayesGLM_EM_old <-
  function(data,
           vertices = NULL,
           faces = NULL,
           mesh = NULL,
           mask = NULL,
           scale_BOLD = TRUE,
           scale_design = TRUE,
           num.threads = 4,
           return_INLA_result = TRUE,
           outfile = NULL,
           verbose = FALSE,
           contrasts = NULL,
           max_iter = 1000,
           tol = 1e-8) {

  # Check to see that the INLA package is installed
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("This function requires the INLA package (see www.r-inla.org/download)")

  #check that only mesh OR vertices+faces supplied
  has_mesh <- !is.null(mesh)
  has_verts_faces <- !is.null(vertices) & !is.null(faces)
  has_howmany <- has_mesh + has_verts_faces
  if(has_howmany != 1) stop('Must supply EITHER mesh OR vertices and faces.')

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data$BOLD)
  n_sess <- length(session_names)

  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  V <- ncol(data$BOLD[[1]]) #number of data locations
  K <- ncol(data$design[[1]]) #number of tasks
  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  }

  if(is.null(mesh)) mesh <- make_mesh(vertices, faces)

  #ID any zero-variance voxels and remove from analysis
  zero_var <- sapply(data, function(x){
    x$BOLD[is.na(x$BOLD)] <- 0 #to detect medial wall locations coded as NA
    x$BOLD[is.nan(x$BOLD)] <- 0 #to detect medial wall locations coded as NaN
    vars <- apply(x$BOLD,2,var)
    return(vars < 1e-6)
  })
  zero_var <- (rowSums(zero_var) > 0) #check whether any vertices have zero variance in any session

  #1. Apply mask to mesh, data and zero_var
  #2. If sum(zero_var) > 0, remove zero_var locations from data and create Amat
  #   Else, let Amat = identity matrix

  if(sum(zero_var) > 0){
    if(!is.null(mask)) mask[zero_var==TRUE] <- 0
    if(is.null(mask)) mask <- !zero_var
  }

  if(!is.null(mask)) {
    mask <- as.logical(mask)
    mesh <- submesh.mesh(mask, mesh)
    mesh$idx$loc <- mesh$idx$loc[!is.na(mesh$idx$loc)]
    for(s in 1:n_sess){
      data[[s]]$BOLD <- data[[s]]$BOLD[,mask]
    }
    V <- sum(mask)
    #zero_var <- zero_var[mask]
  }

  # #remove zero var locations from set of data locations, but leave in the mesh (if no mask provided)
  # Amat <- Diagonal(V, x=1)
  # if(sum(zero_var) > 0){
  #   Amat <- Amat[!zero_var,]
  #   mesh$idx$loc <- mesh$idx$loc[!zero_var]
  #   for(s in 1:n_sess){
  #     data[[s]]$BOLD <- data[[s]]$BOLD[,!zero_var]
  #   }
  # }


  spde <- inla.spde2.matern(mesh)

  #collect data and design matrices
  y_all <- c()
  X_all_list <- NULL
  design <- vector('list', length=n_sess)

  for(s in 1:n_sess){

    #extract and mask BOLD data for current session
    BOLD_s <- data[[s]]$BOLD

    #scale data to represent % signal change (or just center if scale=FALSE)
    BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale_BOLD)
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
  #formula <- make_formula(beta_names = names(betas), repl_names = names(repls), hyper_initial = c(-2,2))
  #formula <- as.formula(formula)

  beta_names <- names(betas)
  repl_names <- names(repls)
  n_beta <- length(names(betas))

  # Hyperparameters
  # hyper_initial <- c(-2,2)
  # hyper_initial <- rep(list(hyper_initial), n_beta)
  # hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

  # formula_vec <- paste0('f(',beta_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
  # formula_vec <- c('y ~ -1', formula_vec)
  # formula_str <- paste(formula_vec, collapse=' + ')
  # formula <- as.formula(formula_str, env = globalenv())

  model_data <- BayesfMRI:::make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)
  od <- organize_data(model_data$y,model_data$X)
  in_binary_mask <- which(mask == 1, arr.ind = T)
  in_binary_mask <- in_binary_mask[,2:1]
  Psi_k <- inla.spde.make.A(mesh,loc = in_binary_mask)
  n <- ncol(Psi_k)

  # Initial values for kappa and tau
  # IMPROVE INITIALS ----
  kappa2 <- tau <- rep(1,K)
  sigma2 <- var(od$y)
  # Q_k <- spde_Q_tau(kappa = kappa,tau = tau,spde = spde)
  Q_k <- spde_Q_phi(kappa2 = kappa2,phi=phi,spde = spde)
  Psi <- Matrix::bdiag(rep(list(Psi_k),K))
  Q_new <- Matrix::bdiag(mapply(spde_Q_phi,kappa2 = kappa2, phi = phi, MoreArgs = list(spde = spde), SIMPLIFY = F))
  Q_new <- Matrix::bdiag(rep(list(Q_k),K))

  -length(model_data$y)/2 * log(var(model_data$y)) - crossprod(model_data$y - mean(model_data$y))/(2*var(model_data$y))

  diff_Theta <- Inf
  step <- 1
  while(norm_diff_Theta > tol) {
    Q <- Q_new
    Sig_inv <- Q + crossprod(model_data$X%*%Psi)
    m <- t(model_data$X%*%Psi)%*%model_data$y / sigma2
    mu <- solve(Sig_inv,m)
    Q_chol <- chol(Q, pivot = T)
    Q_chol_inv <- solve(Q_chol)
    A <- crossprod(model_data$X%*%Psi)
    Sigma <- Q_chol_inv %*% (A/sigma2 + diag(1,length(diag(A)))) %*% t(Q_chol_inv)
    llik <- -length(model_data$y)/2 *  log(sigma2) - crossprod(model_data$y - as.numeric(model_data$X%*%Psi%*%mu))/(2*sigma2)
    cat("Step",step, "kappa =",kappa, "tau =", tau, "sigma^2 =",sigma2 ,"log-likelihood =",llik,"\n")
    for(k in seq(K)) {
      k_inds <- seq(n) + (k-1)*n
      optim_output_k <- optim(
        par = c(kappa[k], tau[k]),
        fn = E_loglik_2,
        Sigma = Sigma[k_inds, k_inds],
        mu = mu[k_inds],
        spde = spde,
        lower = 0,
        method = "L-BFGS-B"
      )
    }
    optim_output <- optim(par = c(kappa,tau), fn = E_loglik_2, Sigma = Sigma, mu = mu, spde = spde, K=K)
    kappa_new <- optim_output$par[1]
    tau_new <- optim_output$par[2]
    Qk_new <- spde_Q(kappa_new,tau_new,spde)
    Q_new <- Matrix::bdiag(rep(list(Qk_new),K))
    sigma2_new <- as.numeric(crossprod(model_data$y) -
      2*crossprod(model_data$y,model_data$X%*%Psi%*%mu) +
      crossprod(model_data$X%*%Psi%*%mu) + sum(colSums(A*Sigma))) /
      length(model_data$y)
    diff_Theta <- c(kappa_new,tau_new, sigma2_new) - c(kappa,tau,sigma2)
    norm_diff_Theta <- sqrt(sum(diff_Theta^2))
    kappa <- kappa_new
    tau <- tau_new
    sigma2 <- sigma2_new
    step <- step+1
  }


  # #estimate model using INLA
  # cat('\n ...... estimating model with INLA')
  # system.time(INLA_result <- BayesfMRI:::estimate_model(formula=formula, data=model_data, A=model_data$X, spde, prec_initial=1, num.threads=num.threads, verbose=verbose, contrasts = contrasts))
  # cat('\n ...... model estimation completed')

  #extract useful stuff from INLA model result
  # beta_estimates <- BayesfMRI:::extract_estimates(object=INLA_result, session_names=session_names, mask=mask) #posterior means of latent task field
  # theta_posteriors <- BayesfMRI:::get_posterior_densities(object=INLA_result, spde, beta_names) #hyperparameter posterior densities

  #extract stuff needed for group analysis
  # mu.theta <- INLA_result$misc$theta.mode
  # Q.theta <- solve(INLA_result$misc$cov.intern)

  #construct object to be returned
  result <- list(INLA_result = NULL,
                 mesh = mesh,
                 mask = mask,
                 design = design,
                 session_names = session_names,
                 beta_names = beta_names,
                 beta_estimates = beta_estimates,
                 theta_posteriors = theta_posteriors,
                 mu.theta = mu.theta, #for joint group model
                 Q.theta = Q.theta, #for joint group model
                 y = y_all, #for joint group model
                 X = X_all_list, #for joint group model
                 #model_data, #temporary
                 #formula, #temporary
                 call = match.call())

  class(result) <- "BayesGLM"

  if(!is.null(outfile)){
    save(result, file=outfile)
  }

  return(result)
}

#' #' Calculate the SPDE covariance
#' #'
#' #' @param kappa A scalar
#' #' @param tau A scalar
#' #' @param spde An \code{inla.spde2} object containing the information about the
#' #'   mesh structure for the SPDE prior
#' #'
#' #' @return The SPDE prior matrix
#' #' @keywords internal
#' spde_Q_tau <- function(kappa,tau, spde) {
#'   Cmat <- spde$param.inla$M0
#'   Gmat <- (spde$param.inla$M1 + Matrix::t(spde$param.inla$M1))/2
#'   GtCinvG <- spde$param.inla$M2
#'   Q <- tau^2 * (kappa^4 * Cmat + 2*kappa^2 * Gmat + GtCinvG)
#'   return(Q)
#' }

#' Calculate the SPDE covariance
#'
#' @param kappa2 A scalar
#' @param phi A scalar
#' @param spde An \code{inla.spde2} object containing the information about the
#'   mesh structure for the SPDE prior
#'
#' @return The SPDE prior matrix
#' @keywords internal
spde_Q_phi <- function(kappa2, phi, spde) {
  Cmat <- spde$param.inla$M0
  Gmat <- (spde$param.inla$M1 + Matrix::t(spde$param.inla$M1))/2
  GtCinvG <- spde$param.inla$M2
  Q <- (kappa2*Cmat + 2*Gmat + GtCinvG/kappa2) / (4*pi*phi)
  return(Q)
}

# E_loglik_2tau <- function(kappa_tau, Sigma, mu, spde) {
#   Q_new <- spde_Q_tau(kappa = kappa_tau[1], tau = kappa_tau[2], spde = spde)
#   # Q_new <- Matrix::bdiag(rep(list(Qk_new),K))
#   trace_QEww <- sum(colSums(Q_new*Sigma)) + crossprod(mu,Q_new)%*%mu
#   log_det_Q <- sum(log(diag(chol(Q,pivot = T))))
#   return(as.numeric(log_det_Q + trace_QEww/2))
# }

# neg_E_loglik_2phi <- function(kappa2, phi, Sigma, mu, spde) {
#   Q_new <- spde_Q_phi(kappa2 = kappa2, phi = phi, spde = spde)
#   trace_QEww <- sum(colSums(Q_new*Sigma)) + crossprod(mu,Q_new)%*%mu
#   log_det_Q <- sum(log(diag(chol(Q_new,pivot = T))))
#   return(-as.numeric(log_det_Q + trace_QEww/2))
# }

#' Gives the portion of the Q matrix independent of phi
#'
#' @param kappa2 scalar
#' @param spdean spde object
#'
#' @return a dgCMatrix
#' @keywords internal
Q_prime <- function(kappa2, spde) {
  Cmat <- spde$param.inla$M0
  Gmat <- (spde$param.inla$M1 + Matrix::t(spde$param.inla$M1))/2
  GtCinvG <- spde$param.inla$M2
  Q <- (kappa2*Cmat + 2*Gmat + GtCinvG/kappa2)
  return(Q)
}

#' The negative of the objective function for kappa
#'
#' @param kappa2 scalar
#' @param spde an spde object
#' @param phi scalar
#' @param Sigma dgCMatrix
#' @param mu dgeMatrix
#'
#' @return a scalar
#' @keywords internal
neg_kappa_fn <- function(kappa2, spde, phi, Sigma, mu) {
  Qp <- Q_prime(kappa2, spde)
  log_det_Q <- sum(log(diag(chol(Qp,pivot = T))))
  trace_QEww <- sum(colSums(Qp*Sigma)) + crossprod(mu,Qp)%*%mu
  out <- (trace_QEww / (4*pi*phi) - log_det_Q)@x
  return(out)
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
#' @inheritParams num.threads_Param
#' @inheritParams return_INLA_result_Param_TRUE
#' @param outfile File name where results will be written (for use by
#'  \code{BayesGLM2}).
#' @inheritParams verbose_Param_inla
#' @inheritParams contrasts_Param_inla
#' @inheritParams avg_betas_over_sessions_Param
#' @param trim_INLA (logical) should the \code{INLA_result} objects within the
#'   result be trimmed to only what is necessary to use the
#'   \code{id_activations} function? Default value is \code{TRUE}.
#'
#' @return A list containing...
#'
#' @importFrom INLA inla.spde2.matern inla.pardiso.check inla.setOption inla.make.lincombs
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#'
#' @export
BayesGLM_EM <- function(
  data, vertices = NULL, faces = NULL, mesh = NULL, mask = NULL,
  scale_BOLD=TRUE, scale_design = TRUE, num.threads=4, return_INLA_result=TRUE,
  outfile = NULL, verbose=FALSE, contrasts = NULL,
  avg_betas_over_sessions = FALSE, trim_INLA = TRUE
  ) {
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

  if(sum(zero_var) > 0){
    if(!is.null(mask)) mask[zero_var==TRUE] <- 0
    if(is.null(mask)) mask <- !zero_var
  }

  if(!is.null(mask)) {
    mask <- as.logical(mask)
    mesh <- excursions::submesh.mesh(mask, mesh)
    mesh$idx$loc <- mesh$idx$loc[!is.na(mesh$idx$loc)]
    for(s in 1:n_sess){
      data[[s]]$BOLD <- data[[s]]$BOLD[,mask]
    }
    V <- sum(mask)
  }

  spde <- inla.spde2.matern(mesh)

  #collect data and design matrices
  y_all <- c()
  X_all_list <- NULL
  design <- vector('list', length=n_sess)

  for(s in 1:n_sess){

    #extract and mask BOLD data for current session
    BOLD_s <- data[[s]]$BOLD

    #scale data to represent % signal change (or just center if scale=FALSE)
    BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale_BOLD)
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
  Psi_k <- inla.spde.make.A(spde)
  Psi <- Matrix::bdiag(rep(list(Psi_k),K))
  A <- crossprod(model_data$X%*%Psi)
  # Initial values for kappa and tau
  kappa2_new <- rep(4,K) # This is a value that matches BayesGLM
  phi_new <- 1 / (4*pi*kappa2_new*4) # This is a value that matches BayesGLM
  sigma2_new <- var(model_data$y)
  Q_k <- mapply(spde_Q_phi, kappa2 = kappa2_new, phi = phi_new, MoreArgs = list(spde = spde), SIMPLIFY = F)
  Q_new <- Matrix::bdiag(Q_k)
  Sig_inv <- Q_new + A/sigma2_new
  m <- t(model_data$X%*%Psi)%*%model_data$y / sigma2_new
  mu <- solve(Sig_inv,m)
  Q_chol <- chol(Q, pivot = T)
  Q_chol_inv <- solve(Q_chol)
  Sigma_new <- Q_chol_inv %*% (A/sigma2 + diag(1,length(diag(A)))) %*% t(Q_chol_inv)
  # mean model log-likelihood (for checking)
  llik_new <- -length(model_data$y)/2 * log(sigma2_new) -
    crossprod(model_data$y - mean(model_data$y))/(2*sigma2_new)
  cat("Log-likelihood for the mean model:", llik_new,"\n")
  llik <- -Inf
  # > EM algorithm ----
  cat(".... performing the EM algorithm\n")
  # diff_Theta <- Inf
  step <- 1
  while(llik_new > llik) {
    Q <- Q_new
    kappa2 <- kappa2_new
    phi <- phi_new
    sigma2 <- sigma2_new
    Sigma <- Sigma_new
    llik <- llik_new
    for(k in seq(K)) {
      k_inds <- seq(V) + (k-1)*V
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
    Sigma <- Q_chol_inv %*% (A/sigma2 + diag(1,length(diag(A)))) %*% t(Q_chol_inv)


    sigma2_new <-
      as.numeric(crossprod(model_data$y) -
                   2*crossprod(model_data$y,model_data$X%*%Psi%*%mu) +
                   crossprod(model_data$X%*%Psi%*%mu) +
                   sum(colSums(A*Sigma))) / length(model_data$y)

    llik_new <- -length(model_data$y)/2 *  log(sigma2) -
      crossprod(model_data$y - as.numeric(model_data$X%*%Psi%*%mu))/(2*sigma2)
    cat("Step",step, "kappa^2 =",kappa2_new, "phi =", phi_new, "sigma^2 =",sigma2_new ,
        "log-likelihood =",llik_new,"\n")
    step <- step+1
  }
  cat(".... EM algorithm complete!")

  Sig_inv <- Q + A/sigma2
  m <- t(model_data$X%*%Psi)%*%model_data$y / sigma2
  mu <- solve(Sig_inv,m)
  Q_chol <- chol(Q, pivot = T)
  Q_chol_inv <- solve(Q_chol)
  Sigma_new <- Q_chol_inv %*% (A/sigma2 + diag(1,length(diag(A)))) %*% t(Q_chol_inv)

  beta_estimates <- matrix(mu@x,nrow = V, ncol = K)
  colnames(beta_estimates) <- beta_names
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
