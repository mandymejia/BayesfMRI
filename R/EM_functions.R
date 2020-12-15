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
#' @param EM_method Either "joint" or "separate" for choosing whether covariates
#'   should share hyperparameter values.
#' @param use_SQUAREM (logical) Should the SQUAREM package be used to speed up
#'   convergence?
#' @param tol A small, positive scalar that determines when iterations should be
#'   terminated.  Default is 1e-3. This is used only if \code{use_SQUAREM} is
#'   \code{TRUE}.
#' @param pct_change_limit The numeric threshold for the percent
#'   (NOT probability) change below which the EM algorithm will stop updates
#'   (default = 1). This is used only if \code{use_SQUAREM} is \code{FALSE}.
#' @param num_cores (optional) allows users to specify the number of cores used
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
#'
#' @export
BayesGLMEM <- function(data,
                       vertices = NULL,
                       faces = NULL,
                       mesh = NULL,
                       mask = NULL,
                       scale_BOLD = TRUE,
                       scale_design = TRUE,
                       EM_method = "joint",
                       use_SQUAREM = TRUE,
                       tol = 1e-3,
                       pct_change_limit = 1,
                       num_cores = 1,
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
  if(!is.null(in_mask)) {
    Psi_k <- INLA::inla.spde.make.A(mesh = mesh, loc = in_mask)
    Psi <- Matrix::bdiag(rep(list(Psi_k),K))
    A <- Matrix::crossprod(model_data$X%*%Psi)
  } else {
    Psi <- Matrix::Diagonal(n = ncol(model_data$X))
    A <- Matrix::crossprod(model_data$X)
  }
  # Initial values for kappa and tau
  if(EM_method == "joint") {
    kappa2 <- 4
  } else {
    kappa2 <- rep(4,K)
  }
  # kappa2 <- ifelse(EM_method == "joint", 4, rep(4,K)) # This is a value that matches BayesGLM
  phi <- 1 / (4*pi*kappa2*4) # This is a value that matches BayesGLM
  sigma2 <- var(model_data$y)
  theta <- c(kappa2, phi, sigma2)
  # > Start EM algorithm ----
  if(EM_method == "joint") {
    # >>  Joint update ----
    if(use_SQUAREM) {
      squareem_output <-
        squarem(
          par = theta,
          fixptfn = GLMEM_fixptjoint,
          objfn = GLMEM_objfn,
          control = list(tol = 1e-3, trace = verbose),
          spde = spde,
          model_data = model_data,
          Psi = Psi,
          K = K,
          A = A
        )
      theta_new <- squareem_output$par
      kappa2_new <- theta_new[1]
      phi_new <- theta_new[2]
      sigma2_new <- theta_new[3]
    } else {
      step <- 1
      max_pct_change <- Inf
      while(max_pct_change > pct_change_limit) {
        theta_new <-
          GLMEM_fixptjoint(
            theta = theta,
            spde = spde,
            model_data = model_data,
            Psi = Psi,
            K = K,
            A = A
          )
        kappa2_new <- theta_new[1]
        phi_new <- theta_new[2]
        sigma2_new <- theta_new[3]
        sigma2_pct_change <- 100*abs((sigma2_new - sigma2) / sigma2)
        phi_pct_change <- 100*abs((phi_new - phi) / phi)
        kappa2_pct_change <- 100*abs((kappa2_new - kappa2) / kappa2)
        max_pct_change <- max(sigma2_pct_change,phi_pct_change,kappa2_pct_change)
        if(verbose) {
          cat("Step",step, "kappa^2 (%change) =",kappa2_new,"(",kappa2_pct_change,")", "phi (%change) =", phi_new, "(", phi_pct_change,")", "sigma^2 (%change) =",sigma2_new,"(",sigma2_pct_change,")","\n")
        }
        kappa2 <- kappa2_new
        phi <- phi_new
        sigma2 <- sigma2_new
        theta <- theta_new
        step <- step+1
      }
    }
  } else {
    cat("Performing separate update. \n")
    # >> Separate update ----
    if(use_SQUAREM) {
      squareem_output <-
        squarem(
          par = theta,
          fixptfn = GLMEM_fixptseparate,
          objfn = GLMEM_objfn,
          control = list(tol = 1e-3, trace = verbose),
          spde = spde,
          model_data = model_data,
          Psi = Psi,
          K = K,
          A = A,
          num_cores = num_cores
        )
      theta_new <- squareem_output$par
      kappa2_new <- theta_new[seq(K)]
      phi_new <- theta_new[seq(K) + K]
      sigma2_new <- theta_new[(2*K + 1)]
    } else {
      step <- 1
      max_pct_change <- Inf
      while(max_pct_change > pct_change_limit) {
        cat("Fixed point function. \n")
        theta_new <-
          GLMEM_fixptseparate(
            theta = theta,
            spde = spde,
            model_data = model_data,
            Psi = Psi,
            K = K,
            A = A,
            num_cores = num_cores
          )
        kappa2_new <- theta_new[seq(K)]
        phi_new <- theta_new[seq(K) + K]
        sigma2_new <- theta_new[(2*K + 1)]
        sigma2_pct_change <- 100*abs((sigma2_new - sigma2) / sigma2)
        phi_pct_change <- 100*abs((phi_new - phi) / phi)
        kappa2_pct_change <- 100*abs((kappa2_new - kappa2) / kappa2)
        max_pct_change <- max(sigma2_pct_change,phi_pct_change,kappa2_pct_change)
        if(verbose) {
          cat("Step",step, "kappa^2 (%change) =",kappa2_new,"(",kappa2_pct_change,")", "phi (%change) =", phi_new, "(", phi_pct_change,")", "sigma^2 (%change) =",sigma2_new,"(",sigma2_pct_change,")","\n")
        }
        kappa2 <- kappa2_new
        phi <- phi_new
        sigma2 <- sigma2_new
        theta <- theta_new
        step <- step+1
      }
    }
  }
  # > End EM algorithm ----
  cat(".... EM algorithm complete!")

  Qk_new <- mapply(spde_Q_phi,kappa2 = kappa2_new, phi = phi_new,
                   MoreArgs = list(spde=spde), SIMPLIFY = F)
  if(EM_method == "joint") {
    Q <- Matrix::bdiag(rep(Qk_new,K))
  } else {
    Q <- Matrix::bdiag(Qk_new)
  }
  Sig_inv <- Q + A/sigma2_new
  m <- Matrix::t(model_data$X%*%Psi)%*%model_data$y / sigma2_new
  mu <- INLA::inla.qsolve(Sig_inv,m)
  Sigma_new <- INLA::inla.qsolve(Sig_inv, Matrix::Diagonal(n = nrow(Sig_inv)), method = "solve")

  beta_estimates <- matrix(mu,nrow = spde$n.spde, ncol = K)
  colnames(beta_estimates) <- beta_names
  beta_estimates <- list(beta_estimates)
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
                 mu.theta = mu.theta, #for joint group model
                 y = y_all, #for joint group model
                 X = X_all_list, #for joint group model
                 call = match.call())
  class(result) <- "BayesGLM"
  return(result)
}

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
#' @param EM_method Either "joint" or "separate" for choosing whether covariates
#'   should share hyperparameter values.
#' @param use_SQUAREM (logical) Should the SQUAREM package be used to speed up
#'   convergence?
#' @param GLM_method Either 'Bayesian' for spatial Bayesian GLM only,
#'   'classical' for the classical GLM only, or 'both' to return both classical
#'   and Bayesian estimates of task activation.
#' @param pct_change_limit The numeric threshold for the percent
#'   (NOT probability) change below which the EM algorithm will stop updates
#'   (default = 1).
#' @param num_cores (optional) allows users to specify the number of cores used
#'   to work in parallel across the different tasks
#' @param session_names (Optional) A vector of names corresponding to each
#'   session.
#' @param outfile (Optional) File name (without extension) of output file for
#'   BayesGLMEM result to use in Bayesian group modeling.
#' @inheritParams verbose_Param_direct_TRUE
#'
#' @importFrom INLA inla.spde2.matern inla.qinv inla.qsolve inla.spde.make.A
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#' @importFrom SQUAREM squarem
#' @importFrom utils head
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
  EM_method = "joint",
  use_SQUAREM = TRUE,
  GLM_method = 'both',
  pct_change_limit = 1,
  num_cores = 1,
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

  # cat(str(binary_mask),"\n")

  ### FIT GLM(s)

  if(do_classical) classicalGLM_out <- classicalGLM(session_data,
                                                    scale_BOLD=scale_BOLD,
                                                    scale_design = scale_design)
  # > Call BayesGLMEM ----
  if(do_Bayesian) BayesGLM_out <- BayesGLMEM(session_data,
                                           mesh = mesh,
                                           mask = binary_mask,
                                           scale_BOLD=scale_BOLD,
                                           scale_design = scale_design,
                                           EM_method = EM_method,
                                           use_SQUAREM = use_SQUAREM,
                                           pct_change_limit = pct_change_limit,
                                           num_cores = num_cores,
                                           outfile = outfile,
                                           verbose=verbose)


  # Extract the point estimates
  # cat(str(convert_mat_A),"\n")
  # cat(str(BayesGLM_out$beta_estimates),"\n")
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
      Bayes_slice[[ss]] <- sapply(seq(num_tasks), function(tn) {
        image_coef <- binary_mask
        image_coef[image_coef == 1] <- point_estimates[[ss]][,tn]
        image_coef[binary_mask == 0] <- NA
        return(image_coef)
      },simplify = F)
      names(Bayes_slice[[ss]]) <- BayesGLM_out$beta_names
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

#' BayesGLM Expectation Maximization for CIFTI
#'
#' Performs spatial Bayesian GLM expectation maximization on the cortical
#' surface for fMRI task activation
#'
#' @section Connectome Workbench Requirement:
#'  This function uses a system wrapper for the 'wb_command' executable. The
#'  user must first download and install the Connectome Workbench, available
#'  from https://www.humanconnectome.org/software/get-connectome-workbench .
#'  The \code{wb_path} argument is the full file path to the Connectome
#'  Workbench folder. (The full file path to the 'wb_cmd' executable also
#'  works.)
#'
# @section Label Levels:
#  \code{xifti$meta$subcort$labels} is a factor with the following levels:
#
#  \enumerate{
#    \item{Cortex-L}
#    \item{Cortex-R}
#    \item{Accumbens-L}
#    \item{Accumbens-R}
#    \item{Amygdala-L}
#    \item{Amygdala-R}
#    \item{Brain Stem}
#    \item{Caudate-L}
#    \item{Caudate-R}
#    \item{Cerebellum-L}
#    \item{Cerebellum-R}
#    \item{Diencephalon-L}
#    \item{Diencephalon-R}
#    \item{Hippocampus-L}
#    \item{Hippocampus-R}
#    \item{Pallidum-L}
#    \item{Pallidum-R}
#    \item{Putamen-L}
#    \item{Putamen-R}
#    \item{Thalamus-L}
#    \item{Thalamus-R}
#  }
#
#  These correspond to the same structures as given by
#  \code{ft_read_cifti} in the \code{cifti-matlab} MATLAB toolbox.
#
#' @param cifti_fname File path (or vector thereof, for multiple-session modeling) of CIFTI-format fMRI timeseries data (*.dtseries.nii).
#' @param surfL_fname File path of GIFTI-format left cortical surface (*.surf.gii). Must be provided if brainstructures includes "left" and GLM_method is "Bayesian" or "both".
#' @param surfR_fname File path of GIFTI-format right cortical surface (*.surf.gii). Must be provided if brainstructures includes "right" and GLM_method is "Bayesian" or "both".
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface) and/or \code{"right"} (right
#'  cortical surface). Default: \code{c("left","right")} (entire cortical surface).
#'  Note that the subcortical models have not yet been implemented.
#' @param wb_path (Optional) Path to Connectome Workbench folder or executable.
#'  If not provided, should be set with
#'  \code{ciftiTools.setOption("wb_path", "path/to/workbench")}.
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
#' @param nuisance (Optional) A TxJ matrix of nuisance signals (or list of such matrices, for multiple-session modeling).
#' @param nuisance_include (Optional) Additional nuisance covariates to include.  Default is 'drift' (linear and quadratic drift terms) and 'dHRF' (temporal derivative of each column of design matrix).
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @param EM_method Either "joint" or "separate" for choosing whether covariates
#'   should share hyperparameter values.
#' @param use_SQUAREM (logical) Should the SQUAREM package be used to speed up
#'   convergence?
#' @param pct_change_limit The numeric threshold for the percent
#'   (NOT probability) change below which the EM algorithm will stop updates
#'   (default = 1).
#' @param num_cores (optional) allows users to specify the number of cores used
#'   to work in parallel across the different tasks
#' @param GLM_method Either 'Bayesian' for spatial Bayesian GLM EM only, 'classical' for the classical GLM only, or 'both' to return both classical and Bayesian estimates of task activation.
#' @param session_names (Optional) A vector of names corresponding to each
#'   session.
#' @param resamp_res The number of vertices to which each cortical surface should be resampled, or NULL if no resampling is to be performed. For computational feasibility, a value of 10000 or lower is recommended.
#' @inheritParams verbose_Param_inla
#' @param outfile (Optional) File name (without extension) of output file for
#'  \code{"BayesGLM"} result to use in Bayesian group modeling.
#'  \code{"_left.rds"} or \code{"_right.rds"} will be appended for the left
#'  cortex and right cortex results, respectively. Default: \code{NULL}
#'  (do not save the results to any file).
#'
#' @return An object of class \code{"BayesGLM"}, a list containing...
#'
#' @importFrom ciftiTools read_cifti resample_gifti as.xifti
#' @importFrom matrixStats rowVars rowSums2 colVars
#' @importFrom INLA inla.spde2.matern inla.qinv inla.qsolve
#' @importFrom excursions submesh.mesh
#' @importFrom SQUAREM squarem
#'
#' @export
BayesGLMEM_cifti <- function(cifti_fname,
                           surfL_fname=NULL, surfR_fname=NULL,
                           brainstructures=c('left','right'),
                           wb_path=NULL,
                           design=NULL, onsets=NULL, TR=NULL,
                           nuisance=NULL, nuisance_include=c('drift','dHRF'),
                           scale_BOLD=TRUE, scale_design=TRUE,
                           EM_method = "joint",
                           use_SQUAREM = TRUE,
                           pct_change_limit = 1,
                           num_cores = 1,
                           GLM_method='both',
                           session_names=NULL,
                           resamp_res=10000,
                           verbose=FALSE,
                           outfile=NULL){

  do_Bayesian <- (GLM_method %in% c('both','Bayesian'))
  do_classical <- (GLM_method %in% c('both','classical'))

  # Check that arguments are compatible
  brainstructures <- ciftiTools:::match_input(
    brainstructures, c("left","right"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }
  do_left <- ('left' %in% brainstructures)
  do_right <- ('right' %in% brainstructures)
  do_sub <- FALSE

  if(do_left & is.null(surfL_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  if(do_right & is.null(surfR_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  if((is.null(design) + is.null(onsets)) != 1) stop('design OR onsets must be provided, but not both')
  if(!is.null(onsets) & is.null(TR)) stop('Please provide TR if onsets provided')

  # Name sessions and check compatibility of multi-session arguments
  n_sess <- length(cifti_fname)
  if(n_sess==1){
    if(is.null(session_names)) session_names <- 'single_session'
    if(!is.null(design)) design <- list(design)
    if(!is.null(onsets)) onsets <- list(onsets)
    if(!is.null(nuisance)) nuisance <- list(nuisance)
  } else {
    if(is.null(session_names)) session_names <- paste0('session', 1:n_sess)
    if(!is.null(design)){ if(length(design) != n_sess) stop('If multiple sessions provided (because cifti_fname is a vector), design must be a list of length equal to the number of sessions (or NULL, if onsets provided).') }
    if(!is.null(onsets)){ if(length(onsets) != n_sess) stop('If multiple sessions provided (because cifti_fname is a vector), onsets must be a list of length equal to the number of sessions (or NULL, if design provided).') }
    if(!is.null(nuisance)){ if(length(nuisance) != n_sess) stop('If multiple sessions provided (because cifti_fname is a vector), nuisance must be a list of length equal to the number of sessions (or NULL).') }
  }
  if(length(session_names) != n_sess) stop('If session_names is provided, it must be of the same length as cifti_fname')

  cat('\n SETTING UP DATA \n')

  ### For each session, separate the CIFTI data into left/right/sub and read in files
  if(do_left) cifti_left <- vector('list', n_sess)
  if(do_right) cifti_right <- vector('list', n_sess)

  if(is.null(design)) {
    make_design <- TRUE
    design <- vector('list', length=n_sess)
  }

  for(ss in 1:n_sess){

    cat(paste0('    Reading in data for session ', ss,'\n'))

    if(ss==1){
      cifti_ss <- read_cifti(
        cifti_fname[ss],
        surfL_fname=surfL_fname, surfR_fname=surfR_fname,
        brainstructures=brainstructures,
        resamp_res=resamp_res,
        wb_path=wb_path
      )
      if(do_left) surf_left <- cifti_ss$surf$cortex_left
      if(do_right) surf_right <- cifti_ss$surf$cortex_right
    } else {
      cifti_ss <- read_cifti(
        cifti_fname[ss],
        brainstructures=brainstructures,
        resamp_res=resamp_res,
        wb_path=wb_path
      )
    }

    if(do_left) {
      cifti_left[[ss]] <- matrix(NA, nrow=length(cifti_ss$meta$cortex$medial_wall_mask$left), ncol=ncol(cifti_ss$data$cortex_left))
      cifti_left[[ss]][cifti_ss$meta$cortex$medial_wall_mask$left,] <- cifti_ss$data$cortex_left
      ntime <- ncol(cifti_left[[ss]])
    }
    if(do_right) {
      cifti_right[[ss]] <- matrix(NA, nrow=length(cifti_ss$meta$cortex$medial_wall_mask$right), ncol=ncol(cifti_ss$data$cortex_right))
      cifti_right[[ss]][cifti_ss$meta$cortex$medial_wall_mask$right,] <- cifti_ss$data$cortex_right
      ntime <- ncol(cifti_right[[ss]])
    }

    if(make_design){
      cat(paste0('    Constructing design matrix for session ', ss, '\n'))
      design[[ss]] <- make_HRFs(onsets[[ss]], TR=TR, duration=ntime)
    }

  }

  ### Check that design matrix names consistent across sessions
  if(n_sess > 1){
    tmp <- sapply(design, colnames)
    tmp <- apply(tmp, 1, function(x) length(unique(x)))
    if(max(tmp) > 1) stop('task names must match across sessions for multi-session modeling')
  }

  cat('\n RUNNING MODELS \n')

  classicalGLM_left <- classicalGLM_right <- classicalGLM_vol <- NULL
  BayesGLM_left <- BayesGLM_right <- BayesGLM_vol <- NULL

  ### FORMAT DESIGN MATRIX
  for(ss in 1:n_sess){
    if(scale_design){
      design[[ss]] <- scale_design_mat(design[[ss]])
    } else {
      design[[ss]] <- scale(design[[ss]], scale=FALSE) #center design matrix to eliminate baseline
    }
  }

  ### ADD ADDITIONAL NUISANCE REGRESSORS
  for(ss in 1:n_sess){
    ntime <- nrow(design[[ss]])
    if('drift' %in% nuisance_include){
      drift <- (1:ntime)/ntime
      if(!is.null(nuisance)) nuisance[[ss]] <- cbind(nuisance[[ss]], drift, drift^2) else nuisance[[ss]] <- cbind(drift, drift^2)
    }
    if('dHRF' %in% nuisance_include){
      dHRF <- gradient(design[[ss]])
      if(!is.null(nuisance)) nuisance[[ss]] <- cbind(nuisance[[ss]], dHRF) else nuisance[[ss]] <- dHRF
    }
  }

  scale_design <- FALSE # This is done to prevent double-scaling in the BayesGLM function

  ### LEFT HEMISPHERE
  if(do_left){

    cat('\n ... LEFT CORTEX \n')

    #set up mesh
    verts_left <- surf_left$vertices #first surface is used for modeling
    faces_left <- surf_left$faces

    #set up session list
    session_data <- vector('list', n_sess)
    names(session_data) <- session_names
    for(ss in 1:n_sess){
      sess <- list(BOLD = t(cifti_left[[ss]]), design=design[[ss]])
      if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
      session_data[[ss]] <- sess
    }

    ### FIT GLM(s)

    if(!is.null(outfile)) {
      if (endsWith(outfile, ".rds")) {
        outfile_left <- gsub(".rds$", "_left.rds", outfile)
      } else {
        outfile_left <- paste0(outfile, "_left.rds")
      }
    } else {
      outfile_left <- NULL
    }

    if(do_classical) classicalGLM_left <- classicalGLM(session_data,
                                                       scale_BOLD=scale_BOLD,
                                                       scale_design = scale_design)
    # > Call BayesGLMEM (left) ----
    if(do_Bayesian) BayesGLM_left <- BayesGLMEM(session_data,
                                              vertices = verts_left,
                                              faces = faces_left,
                                              scale_BOLD=scale_BOLD,
                                              scale_design = scale_design,
                                              EM_method = EM_method,
                                              use_SQUAREM = use_SQUAREM,
                                              pct_change_limit = pct_change_limit,
                                              num_cores = num_cores,
                                              outfile = outfile,
                                              verbose=verbose)

  }


  ### RIGHT HEMISPHERE
  if(do_right){

    cat('\n ... RIGHT CORTEX \n')

    #set up mesh
    #surf_right <- readGIfTI(surfR_fname)$data
    verts_right <- surf_right$vertices #first surface is used for modeling
    faces_right <- surf_right$faces
    #if(min(faces_right)==0) faces_right <- faces_right + 1

    #set up session list
    session_data <- vector('list', n_sess)
    names(session_data) <- session_names
    for(ss in 1:n_sess){
      sess <- list(BOLD = t(cifti_right[[ss]]), design=design[[ss]])
      if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
      session_data[[ss]] <- sess
    }

    ### FIT GLM

    if(!is.null(outfile)) {
      if (endsWith(outfile, ".rds")) {
        outfile_right <- gsub(".rds$", "_right.rds", outfile)
      } else {
        outfile_right <- paste0(outfile, "_right.rds")
      }
    } else {
      outfile_right <- NULL
    }

    if(do_classical) classicalGLM_right <- classicalGLM(session_data,
                                                        scale_BOLD=scale_BOLD,
                                                        scale_design = scale_design)
    # > Call BayesGLMEM (right) ----
    if(do_Bayesian) BayesGLM_right <- BayesGLMEM(session_data,
                                                vertices = verts_right,
                                                faces = faces_right,
                                                scale_BOLD=scale_BOLD,
                                                scale_design = scale_design,
                                                EM_method = EM_method,
                                                use_SQUAREM = use_SQUAREM,
                                                pct_change_limit = pct_change_limit,
                                                num_cores = num_cores,
                                                outfile = outfile,
                                                verbose=verbose)
  }

  ### CONSTRUCT BETA ESTIMATES AS CIFTI OBJECTS

  cat('\n PUTTING RESULTS IN CIFTI FORMAT \n')

  classicalGLM_cifti <- BayesGLM_cifti <- vector('list', n_sess)
  names(classicalGLM_cifti) <- names(BayesGLM_cifti) <- session_names
  for(ss in 1:n_sess){
    if(do_classical){
      classicalGLM_cifti[[ss]] <- as.xifti(
        cortexL = classicalGLM_left[[ss]],
        cortexR = classicalGLM_right[[ss]]
      )
    }
    if(do_Bayesian){
      betas_left <- betas_right <- NULL
      if(do_left) {
        betas_left <- matrix(NA,length(mask),ncol(BayesGLM_left$beta_estimates[[ss]]))
        betas_left[mask,] <- BayesGLM_left$beta_estimates[[ss]]
      }
      if(do_right) {
        betas_right <- matrix(NA,length(mask),ncol(BayesGLM_right$beta_estimates[[ss]]))
        betas_right[mask,] <- BayesGLM_left$beta_estimates[[ss]]
      }

      BayesGLM_cifti[[ss]] <- as.xifti(
        cortexL = betas_left,
        cortexR = betas_right
      )
    }
  }

  if (do_Bayesian) {
    if (do_left) {
      beta_names <- BayesGLM_left$beta_names
    } else {
      beta_names <- BayesGLM_right$beta_names
    }
  } else {
    beta_names <- NULL
  }

  result <- list(session_names = session_names,
                 beta_names = beta_names,
                 betas_Bayesian = BayesGLM_cifti,
                 betas_classical = classicalGLM_cifti,
                 GLMs_Bayesian = list(cortexL = BayesGLM_left,
                                      cortexR = BayesGLM_right),
                 GLMs_classical = list(cortexL = classicalGLM_left,
                                       cortexR = classicalGLM_right),
                 design = design)

  cat('\n DONE! \n')

  class(result) <- "BayesGLM_cifti"
  return(result)
}
