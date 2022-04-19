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
#' @inheritParams avg_sessions_Param
#'
#' @return A list containing...
#'
#' @import MatrixModels
#' @importFrom INLA inla.spde2.matern inla.qsolve
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#' @importFrom methods as
#' @importFrom SQUAREM squarem
#' @importFrom parallel makeCluster parApply
#' @importFrom utils tail
#' @importFrom stats var
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
                       use_SQUAREM = TRUE,
                       tol = NULL,
                       num.threads = 1,
                       outfile = NULL,
                       verbose = FALSE,
                       avg_sessions = TRUE) {

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

  if(is.null(mesh)) mesh <- make_mesh(vertices, faces,use_INLA = FALSE)
  # if(is.null(mesh)) mesh_new <- galerkin_db(faces, vertices, surface = TRUE)

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
      # mesh <- excursions::submesh.mesh(mask, mesh)
      mesh <- submesh(mask, mesh)
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
  # spde <- INLA::inla.spde2.matern(mesh)
  spde <- create_spde_surf(mesh)
  spde <- spde$spde

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
  kappa2 <- 4
  phi <- 1 / (4*pi*kappa2*4) # This is a value that matches BayesGLM
  # Using values based on the classical GLM
  cat("... FINDING BEST GUESS INITIAL VALUES\n")
  beta_hat <- MatrixModels:::lm.fit.sparse(model_data$X, model_data$y)
  res_y <- (model_data$y - model_data$X %*% beta_hat)@x
  sigma2 <- stats::var(res_y)
  beta_hat <- matrix(beta_hat, ncol = K*n_sess)
  rcpp_spde <- create_listRcpp(spde)
  spde <- rcpp_spde
  if(n_sess > 1) {
    task_cols <- sapply(seq(n_sess), function(j) seq(K) + K *(j - 1))
    beta_hat <- apply(task_cols,1,function(x) beta_hat[,x])
  }
  # kappa2_phi_rcpp <-
  #   apply(
  #     beta_hat,
  #     2,
  #     initialKP,
  #     theta = c(kappa2, phi),
  #     spde = rcpp_spde,
  #     n_sess = n_sess,
  #     tol = tol,
  #     verbose = verbose
  #   )
  n_threads <- parallel::detectCores()
  n_threads <- min(n_threads,K,num.threads)
  cl <- parallel::makeCluster(n_threads)
  kappa2_phi_rcpp <-
    parallel::parApply(
      cl = cl,
      beta_hat,
      2,
      initialKP,
      theta = c(kappa2, phi),
      spde = rcpp_spde,
      n_sess = n_sess,
      tol = tol,
      verbose = FALSE
    )
  parallel::stopCluster(cl)
  cat("...... DONE!\n")
  theta <- c(t(kappa2_phi_rcpp), sigma2)
  theta_init <- theta
  Ns <- 50
  cat("... STARTING EM ALGORITHM\n")
  em_output <-
    findTheta(
      theta = theta,
      spde = rcpp_spde,
      y = model_data$y,
      X = model_data$X,
      QK = make_Q(theta, rcpp_spde, n_sess),
      Psi = as(Psi, "dgCMatrix"),
      A = as(A, "dgCMatrix"),
      Ns = 50,
      tol = tol,
      verbose = verbose
    )
  # > End EM algorithm ----
  cat(".... EM algorithm complete!\n")
  kappa2_new <- phi_new <- sigma2_new <- mu <- NULL
  list2env(em_output, envir = environment())
  Qk_new <- mapply(spde_Q_phi,kappa2 = kappa2_new, phi = phi_new,
                   MoreArgs = list(spde=rcpp_spde), SIMPLIFY = F)
  Q <- Matrix::bdiag(Qk_new)
  if(n_sess > 1) Q <- Matrix::bdiag(lapply(seq(n_sess), function(x) Q))
  Sig_inv <- Q + A/sigma2_new
  m <- Matrix::t(model_data$X%*%Psi)%*%model_data$y / sigma2_new
  mu <- Matrix::solve(Sig_inv, m)

  beta_estimates <- matrix(mu,nrow = V_nm, ncol = K*n_sess)
  colnames(beta_estimates) <- rep(beta_names, n_sess)
  beta_estimates <- lapply(seq(n_sess), function(ns) beta_estimates[,(seq(K) + K * (ns - 1))])
  names(beta_estimates) <- session_names
  avg_beta_estimates <- NULL
  if(avg_sessions) avg_beta_estimates <- Reduce(`+`,beta_estimates) / n_sess
  theta_estimates <- c(sigma2_new,c(phi_new,kappa2_new))
  names(theta_estimates) <- c("sigma2",paste0("phi_",seq(K)),paste0("kappa2_",seq(K)))
  #extract stuff needed for group analysis
  tau2_init <- 1 / (4*pi*theta_init[seq(K)]*theta_init[(seq(K) + K)])
  mu.theta_init <- c(log(tail(theta_init,1)), c(rbind(log(sqrt(tau2_init)),log(sqrt(theta_init[seq(K)])))))
  tau2 <- 1 / (4*pi*kappa2_new*phi_new)
  mu.theta <- c(log(1/sigma2_new),c(rbind(log(sqrt(tau2)),log(sqrt(kappa2_new)))))

  # Construct object to be returned
  result <- list(mesh = mesh,
                 mask = mask,
                 design = design,
                 session_names = session_names,
                 beta_names = beta_names,
                 beta_estimates = beta_estimates,
                 avg_beta_estimates = avg_beta_estimates,
                 theta_estimates = theta_estimates,
                 posterior_Sig_inv = Sig_inv, # For excursions
                 mu.theta = mu.theta, #for joint group model
                 mu.theta_init = mu.theta_init, # To examine convergence
                 y = y_all, #for joint group model
                 X = X_all_list, #for joint group model
                 call = match.call())
  class(result) <- "BayesGLM"
  return(result)
}

#' BayesGLM for 3D volume
#'
#' Applies spatial Bayesian GLM to task fMRI data for 3D subcortical volumes
#'
#' The subcortical data is separated into regions, whose sizes range
#'  from approximately 100 voxels to approximately 9000 voxels.  Smaller regions
#'  are grouped together to improve model fit.
#'  The \code{groups_df} argument specifies which regions are grouped together.
#'  This argument should be a data frame with R rows (the number of regions) and
#'  three columns: label, region, and group.
#'  The label column is the numerical identifier of each region; the region
#'  column contains the region names, and the group column contains the model
#'  group assignments (e.g. 1,2,3). Regions to be excluded
#'  from analysis are indicated by NA in the group assignment.
#'
#' @param data A list of sessions, where each session is a list with elements
#'  BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#'  List element names represent session names.
#' @param beta_names (Optional) Names of tasks represented in design matrix
#' @param locations Vx3 matrix of x,y,z coordinates of each voxel
#' @param labels Vector of length V of region labels
#' @param use_SQUAREM (logical) Should the SQUAREM package be used to speed up
#'   convergence?
#' @param ar_order (numeric) Controls prewhitening. If greater than zero, this
#'  should be a number indicating the order of the autoregressive model to use
#'  for prewhitening. If zero, do not prewhiten. Default: \code{0}.
#' @param ar_smooth FWHM parameter for smoothing. Remember that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} or \code{NULL}
#'  to not do any smoothing. Default: \code{0}.
#' @param cifti_data A \code{xifti} object with data and metadata pertaining to
#'   the subcortex. Used to smooth AR coefficients during prewhitening.
#' @param tol If use_SQUAREM == TRUE, an absolute change limit for
#'   when the EM algorithm should be stopped (Default = 1e-3). If use_SQUAREM ==
#'   FALSE, a percent change limit for when the EM algorithm should be stopped
#'   (Default = 1). A value of NULL will result in the default value being used.
#' @param groups_df Data frame indicating the name and model group of each region.  See Details.
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @param outfile File name where results will be written (for use by \code{BayesGLM_grp}).
#' @inheritParams num.threads_Param
#' @inheritParams verbose_Param_inla
#' @inheritParams avg_sessions_Param
#'
#' @return A list containing...
#'
#' @export
BayesGLMEM_vol3D <-
  function(data,
           beta_names = NULL,
           locations,
           labels,
           use_SQUAREM = TRUE,
           ar_order = 0,
           ar_smooth = 0,
           cifti_data = NULL,
           tol = NULL,
           groups_df = NULL,
           scale_BOLD = TRUE,
           scale_design = TRUE,
           outfile = NULL,
           num.threads = 6,
           verbose = FALSE,
           avg_sessions = TRUE) {
  if(is.null(groups_df)) {
    regions <- c('Accumbens-l','Accumbens-r', #3,4 -- 200 voxels - BASAL GANGLIA --> MODEL 1
                 'Amygdala-l','Amygdala-r',   #5,6 -- 600 voxels  IMPORTANT --> MODEL 2
                 'Brain Stem',                #7 -- 3,402 voxels  EXCLUDE
                 'Caudate-L','Caudate-R',     #8,9 -- 1400 voxels  BASAL GANGLIA --> MODEL 1
                 'Cerebellum-L','Cerebellum-R', #10,11 -- 9000 voxels each --> MODELS 3 AND 4
                 'Diencephalon-L','Diencephalon-R', #12,13 -- 1400 voxels  --> MODEL 2
                 'Hippocampus-L','Hippocampus-R', #14,15 -- 1500 voxels  IMPORTANT --> MODEL 2
                 'Pallidum-L','Pallidum-R', #16,17 -- 500 voxels  BASAL GANGLIA --> MODEL 1
                 'Putamen-L','Putamen-R', #18,19 -- 2000 voxels  BASAL GANGLIA --> MODEL 1
                 'Thalamus-L','Thalamus-R' #20,21 -- 2500 voxels IMPORTANT --> MODEL 2
    )
    groups <- c(1,1,2,2,NA,1,1,3,4,2,2,2,2,1,1,1,1,2,2)
    groups_df <- data.frame(label=3:21, region=regions, group=groups)
  }

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
  if(!is_pw) {
    K <- ncol(data[[1]]$design) #number of tasks
  } else {
    K <- ncol(data[[1]]$design) / sum(!is_missing)
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

  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(!is_pw) {
      if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
    } else {
      if(ncol(data[[s]]$design) / V_nm != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
    }
  }

  if(is.null(outfile)){
    warning('No value supplied for outfile, which is required for group modeling (see `help(BayesGLM2)`).')
  }

  if(nrow(locations) != V) stop('The locations argument should have V rows, the number of data locations in BOLD.')
  if(ncol(locations) != 3) stop('The locations argument should have 3 columns indicating the x,y,z coordinate of each data location.')

  # Check prewhitening arguments.
  # if(is.null(ar_order)) ar_order <- 0
  # ar_order <- as.numeric(ar_order)
  # prewhiten <- (ar_order > 0)

  ### Run SPDE object for each group of regions

  regions_set <- unique(labels)
  num_regions <- length(regions_set)
  if(class(groups_df) != 'data.frame') stop('The groups_df argument should be a data frame with the following columns: label, region, group.')
  if(!all.equal(sort(names(groups_df)), sort(c('label','region','group')))) stop('The groups_df argument should be a data frame with the following columns: label, region, group.')
  if(nrow(groups_df) != num_regions) stop('The groups_df argument should contain one row for each region represented in the labels vector.')


  # > Fit model for each group of regions ----

  groups_df <- groups_df[!is.na(groups_df$group),]
  group_set <- unique(groups_df$group)

  beta_estimates_all <- matrix(NA, nrow=V, ncol=K)
  beta_estimates_all <- rep(list(beta_estimates_all), n_sess)
  avg_beta_estimates_all <- matrix(NA, nrow=V, ncol=K)
  theta_estimates_all <- matrix(NA, nrow = length(group_set), ncol = 2*K + 1)
  EM_result_all <- vector('list', length=length(group_set))
  names(EM_result_all) <- paste0('model_group_',group_set)
  spde_all <- vector('list', length=length(group_set))
  names(spde_all) <- paste0('model_group_',group_set)
  # pw_data <- vector("list", length = length(group_set))

  # Prewhiten data
  # if(prewhiten) {
  #   pw_data <-
  #     prewhiten_cifti(
  #       data = data,
  #       mask = NULL,
  #       scale_BOLD = TRUE,
  #       scale_design = TRUE,
  #       ar_order = ar_order,
  #       ar_smooth = ar_smooth,
  #       cifti_data = cifti_data,
  #       brainstructure = "subcortical",
  #       num.threads = num.threads
  #     )
  #   session_data <- pw_data$data
  #   scale_BOLD_sub <- FALSE # Done above
  # }

  for(grp in group_set){
    label_set_grp <- groups_df$label[groups_df$group==grp]
    name_set_grp <- groups_df$region[groups_df$group==grp]
    inds_grp <- as.numeric(labels) %in% label_set_grp
    locs_grp <- locations[inds_grp,]
    labels_grp <- labels[inds_grp]

    cat(paste0('Estimating Model ',grp,' (', paste(name_set_grp, collapse = ', '), ')'),"\n")

    spde_grp <- create_spde_vol3D(locs=locs_grp, labs=labels_grp, lab_set=label_set_grp)
    spde <- spde_grp$spde

    is_pw <- FALSE
    if(ar_order > 0) {
      pw_data <-
        prewhiten_cifti(
          data = data,
          mask = inds_grp,
          scale_BOLD = scale_BOLD,
          scale_design = scale_design,
          ar_order = ar_order,
          ar_smooth = ar_smooth,
          cifti_data = cifti_data,
          brainstructure = "subcortical",
          num.threads = num.threads
        )
      scale_BOLD <- scale_design <- FALSE # These have beend done
      is_pw <- TRUE
      data <- pw_data$data
    }

    #collect data and design matrices
    y_all <- c()
    X_all_list <- NULL

    for(s in 1:n_sess){

      #extract and mask BOLD data for current session
      BOLD_s <- data[[s]]$BOLD[,inds_grp]

      #scale data to represent % signal change (or just center if scale=FALSE)
      if(scale_BOLD) {
        BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale_BOLD, transpose = FALSE)
      }
      if(!is_pw)
        design_s <- scale(data[[s]]$design, scale=scale_design) #center design matrix to eliminate baseline

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
      if(!is_pw) {
        data_org <- organize_data(y_reg, X_reg, transpose = FALSE)
      }
      if(is_pw) {
        data_org <- organize_data_pw(y_reg, X_reg, transpose = FALSE)
      }
      y_vec <- data_org$y
      X_list <- list(data_org$X)
      names(X_list) <- session_names[s]

      y_all <- c(y_all, y_vec)
      X_all_list <- c(X_all_list, X_list)
    }

    #construct betas and repls objects
    replicates_list <- organize_replicates(n_sess=n_sess, beta_names = beta_names, mesh = spde_grp)
    betas <- replicates_list$betas
    repls <- replicates_list$repls

    model_data <- make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)

    # >> Model setup and initials ----
    if(is.null(tol)){
      if(use_SQUAREM) tol <- 1e-3
      if(!use_SQUAREM) tol <- 1
    }
    Psi_k <- spde_grp$Amat
    Psi <- Matrix::bdiag(rep(list(Psi_k),K))
    if(n_sess > 1) Psi <- Matrix::bdiag(rep(list(Psi),n_sess))
    A <- Matrix::crossprod(model_data$X%*%Psi)
    # Initial values for kappa and tau
    # Using values matching BayesGLM
    kappa2 <- 4
    phi <- 1 / (4*pi*kappa2*4) # This is a value that matches BayesGLM
    # Using values based on the classical GLM
    cat("... FINDING BEST GUESS INITIAL VALUES\n")
    beta_hat <- MatrixModels:::lm.fit.sparse(model_data$X, model_data$y)
    res_y <- (model_data$y - model_data$X %*% beta_hat)@x
    sigma2 <- stats::var(res_y)
    beta_hat <- (beta_hat %*% Psi)@x
    beta_hat <- matrix(beta_hat, ncol = K*n_sess)
    rcpp_spde <- create_listRcpp(spde)
    spde <- rcpp_spde
    if(n_sess > 1) {
      task_cols <- sapply(seq(n_sess), function(j) seq(K) + K *(j - 1))
      beta_hat <- apply(task_cols,1,function(x) beta_hat[,x])
    }
    n_threads <- parallel::detectCores()
    n_threads <- min(n_threads,K,num.threads)
    cl <- parallel::makeCluster(n_threads)
    kappa2_phi_rcpp <-
      parallel::parApply(
        cl = cl,
        beta_hat,
        2,
        initialKP,
        theta = c(kappa2, phi),
        spde = rcpp_spde,
        n_sess = n_sess,
        tol = tol,
        verbose = verbose
      )
    parallel::stopCluster(cl)
    cat("...... DONE!\n")
    theta <- c(t(kappa2_phi_rcpp), sigma2)
    theta_init <- theta
    Ns <- 50
    cat("... STARTING EM ALGORITHM\n")
    em_output <-
      findTheta(
        theta = theta,
        spde = rcpp_spde,
        y = model_data$y,
        X = model_data$X,
        QK = make_Q(theta, rcpp_spde, n_sess),
        Psi = as(Psi, "dgCMatrix"),
        A = as(A, "dgCMatrix"),
        Ns = 50,
        tol = tol,
        verbose = verbose
      )
    # > End EM algorithm ----
    cat(".... EM algorithm complete!\n")
    kappa2_new <- phi_new <- sigma2_new <- NULL
    list2env(em_output, envir = environment())
    Qk_new <- mapply(spde_Q_phi,kappa2 = kappa2_new, phi = phi_new,
                     MoreArgs = list(spde=spde), SIMPLIFY = F)
    Q <- Matrix::bdiag(Qk_new)
    if(n_sess > 1) Q <- Matrix::bdiag(lapply(seq(n_sess), function(x) Q))
    Sig_inv <- Q + A/sigma2_new
    m <- Matrix::t(model_data$X%*%Psi)%*%model_data$y / sigma2_new
    mu <- Matrix::solve(Sig_inv, m)
    beta_estimates <- matrix(mu,nrow = nrow(spde$Cmat), ncol = K*n_sess)
    colnames(beta_estimates) <- rep(beta_names, n_sess)
    beta_estimates <- as.matrix(Psi_k %*% beta_estimates)
    beta_estimates <- lapply(seq(n_sess), function(ns) beta_estimates[,(seq(K) + K * (ns - 1))])
    names(beta_estimates) <- session_names
    avg_beta_estimates <- NULL
    if(avg_sessions) avg_beta_estimates <- Reduce(`+`,beta_estimates) / n_sess
    theta_estimates <- c(sigma2_new,c(phi_new,kappa2_new))
    names(theta_estimates) <- c("sigma2",paste0("phi_",seq(K)),paste0("kappa2_",seq(K)))
    #extract stuff needed for group analysis
    tau2_init <- 1 / (4*pi*theta_init[seq(K)]*theta_init[(seq(K) + K)])
    mu.theta_init <- c(log(tail(theta_init,1)), c(rbind(log(sqrt(tau2_init)),log(sqrt(theta_init[seq(K)])))))
    tau2 <- 1 / (4*pi*kappa2_new*phi_new)
    mu.theta <- c(log(1/sigma2_new),c(rbind(log(sqrt(tau2)),log(sqrt(kappa2_new))))) # This is a guess about the order and might be wrong

    #extract beta estimates and project back to data locations for current group
    for(s in 1:n_sess){
      beta_estimates_all[[s]][inds_grp,] <- beta_estimates[[s]]
    }
    if(avg_sessions)
      avg_beta_estimates_all[inds_grp,] <- avg_beta_estimates
    #extract theta estimates and project back to data locations for current group
    theta_estimates_all[grp,] <- mu.theta

    # return EM result for each group
    EM_result_all[[grp]] <- list(mesh = spde_grp,
                   session_names = session_names,
                   beta_names = beta_names,
                   beta_estimates = beta_estimates,
                   avg_beta_estimates = avg_beta_estimates,
                   theta_estimates = theta_estimates,
                   posterior_mu = mu, # For excursions
                   posterior_Sig_inv = Sig_inv, # For excursions
                   mu.theta = mu.theta, #for joint group model
                   mu.theta_init = mu.theta_init, # To examine convergence
                   y = y_all, #for joint group model
                   X = X_all_list, #for joint group model
                   call = match.call())

    # Return the SPDE object for each group
    spde_all[[grp]] <- spde_grp
  }
  # Bring it all together
  result <- list(spde_obj = spde_all,
                 mesh = sapply(spde_all, function(x) x$spde$n.spde, simplify = F),
                 session_names = session_names,
                 EM_result_all = EM_result_all,
                 beta_names = beta_names,
                 beta_estimates = beta_estimates_all,
                 avg_beta_estimates = avg_beta_estimates_all,
                 theta_estimates = theta_estimates_all,
                 call = match.call()
  )
  class(result) <- "BayesGLM"
  return(result)
}

