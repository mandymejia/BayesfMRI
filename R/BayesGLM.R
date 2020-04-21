#' Performs whole-brain spatial Bayesian GLM for fMRI task activation
#'
#' @param fname_cifti File path (or vector thereof, for multiple-session modeling) of CIFTI-format fMRI timeseries data (*.dtseries.nii).
#' @param fname_gifti_left File path of GIFTI-format left cortical surface (*.surf.gii). Must be provided if brainstructures includes "left" and GLM_method is "Bayesian" or "both".
#' @param fname_gifti_right File path of GIFTI-format right cortical surface (*.surf.gii). Must be provided if brainstructures includes "right" and GLM_method is "Bayesian" or "both".
#' @param brainstructures A vector indicating which brain structure(s) to model: 'left' (left cortical surface), 'right' (right cortical surface), and/or 'subcortical' (subcortical and cerebellar gray matter)
#' @param vol_regions A vector indicating which subcortical brain regions (3-21) to model. Default is to exclude brainstem (region 7).
#' @param wb_cmd Path to Connectome Workbench executable file, ending in 'wb_command' (Mac/linux) or 'wb_command.exe' (Windows).
#' @param design A TxK task design matrix (or list of such matrices, for multiple-session modeling) with column names representing tasks. Each column represents the expected BOLD response due to each task, a convolution of the hemodynamic response function (HRF) and the task stimulus.  Must be provided if and only if onsets=NULL.
#' @param onsets A matrix of onsets and durations for each task (or a list of such matrices, for multiple-session modeling).  Must be provided if and only if design=NULL.
#' @param nuisance (Optional) A TxJ matrix of nuisance signals (or list of such matrices, for multiple-session modeling).
#' @param GLM_method Either 'Bayesian' for spatial Bayesian GLM only, 'classical' for the classical GLM only, or 'both' to return both classical and Bayesian estimates of task activation.
#' @param session_names (Optional) A vector of names corresponding to each session.
#'
#' @return An object of class BayesGLM, a list containing ...
#' @export
#' @importFrom ciftiTools cifti_read_separate
#' @importFrom matrixStats rowVars rowSums2
#'
#' @details This function uses a system wrapper for the 'wb_command' executable. The user must first download and install the Connectome Workbench,
#' available from https://www.humanconnectome.org/software/get-connectome-workbench. The 'wb_cmd' argument is the full file path to the 'wb_command' executable file.
#'
#' The subcortical brain structure labels range from 3-21 and represent:
#' 3 Accumbens-L
#' 4 Accumbens-R
#' 5 Amygdala-L
#' 6 Amygdala-R
#' 7 Brain Stem
#' 8 Caudate-L
#' 9 Caudate-R
#' 10 Cerebellum-L
#' 11 Cerebellum-R
#' 12 Diencephalon-L
#' 13 Diencephalon-R
#' 14 Hippocampus-L
#' 15 Hippocampus-R
#' 16 Pallidum-L
#' 17 Pallidum-R
#' 18 Putamen-L
#' 19 Putamen-R
#' 20 Thalamus-L
#' 21 Thalamus-R
#'
BayesGLM <- function(fname_cifti, fname_gifti_left=NULL, fname_gifti_right=NULL, brainstructures=c('left','right','subcortical'), vol_regions=c(3:6,8:21), wb_cmd, design=NULL, onsets=NULL, nuisance=NULL, GLM_method='both', session_names=NULL){

  # Check that arguments are compatible
  do_left <- ('left' %in% brainstructures)
  do_right <- ('right' %in% brainstructures)
  do_sub <- ('subcortical' %in% brainstructures)
  if(do_left & is.null(fname_gifti_left)) stop('fname_gifti_left must be provided if brainstructures includes "left"')
  if(do_right & is.null(fname_gifti_right)) stop('fname_gifti_left must be provided if brainstructures includes "left"')
  if((is.null(design) + is.null(onsets)) != 1) stop('design OR onsets must be provided, but not both')
  if(do_sub){
    if(length(unique(vol_regions)) != length(vol_regions)) stop('vol_regions must contain no repeated values.')
    if(min(is.element(vol_regions, 3:21))==0) stop('vol_regions must include only integer values between 3 and 21.')
  }

  # Name sessions and check compatibility of multi-session arguments
  n_sess <- length(fname_cifti)
  if(n_sess==1){
    if(is.null(session_names)) session_names <- 'single_session'
    fname_cifti <- list(fname_cifti)
    if(!is.null(design)) design <- list(design)
    #if(!is.null(onsets)) onsets <- list(onsets)
    if(!is.null(nuisance)) nuisance <- list(nuisance)
  } else {
    if(is.null(session_names)) session_names <- paste0('session', 1:n_sess)
    if(!is.null(design)){ if(length(design) != n_sess) stop('If multiple sessions provided (because fname_cifti is a vector), design must be a list of length equal to the number of sessions (or NULL, if onsets provided).') }
    #if(!is.null(onsets)){ if(length(onsets) != n_sess) stop('If multiple sessions provided (because fname_cifti is a vector), onsets must be a list of length equal to the number of sessions (or NULL, if design provided).') }
    if(!is.null(nuisance)){ if(length(nuisance) != n_sess) stop('If multiple sessions provided (because fname_cifti is a vector), nuisance must be a list of length equal to the number of sessions (or NULL).') }
  }
  if(length(session_names) != n_sess) stop('If session_names is provided, it must be of the same length as fname_cifti')


  ### For each session, separate the CIFTI data into left/right/sub and read in files
  if(do_left) cifti_left <- vector('list', n_sess)
  if(do_right) cifti_right <- vector('list', n_sess)
  if(do_sub) nifti_data <- nifti_labels <- vector('list', n_sess)
  for(ss in 1:n_sess){
    cifti_ss <- cifti_read_separate(fname_cifti[[ss]], brainstructures=brainstructures)
    if(do_left) cifti_left[[ss]] <- cifti_ss$CORTEX_LEFT
    if(do_right) cifti_right[[ss]] <- cifti_ss$CORTEX_RIGHT
    if(do_sub) nifti_data[[ss]] <- cifti_ss$VOL
    if(do_sub & ss==1) nifti_labels[[ss]] <- cifti_ss$LABELS
  }
  #check that labels are the same across all sessions
  if(do_sub) {
    if(n_sess > 1) {
      tmp <- sapply(nifti_labels, function(x) {all.equal(x,nifti_labels[[1]])})
      if(min(tmp)==0) stop('Subcortical labels must match across all sessions in cifti data. Check compatibility of cifti files.')
    }
    nifti_labels <- nifti_labels[[1]]
  }


  ### TO DO: DEAL WITH MEDIAL WALL LOCATIONS.  BAYESGLM WILL NOT WORK WITH NAS.

  classicalGLM_left <- classicalGLM_right <- classicalGLM_vol <- NULL
  BayesGLM_left <- BayesGLM_right <- BayesGLM_vol <- NULL

  ### LEFT HEMISPHERE
  if(do_left){

    #set up mesh
    surf_left <- readGIfTI(fname_gifti_left)$data
    verts_left <- surf_left$pointset
    faces_left <- surf_left$triangle + 1

    #set up session list
    session_data <- vector('list', n_sess)
    names(session_data) <- session_names
    for(ss in 1:n_sess){
       sess <- list(BOLD = t(cifti_left[[ss]]), design=design[[ss]])
      if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
      session_data[[ss]] <- sess
    }

    ### FIT GLM

    if(GLM_method %in% c('classical','both')) classicalGLM_left <- classicalGLM(session_data) else classicalGLM_left <- NULL
    if(GLM_method %in% c('Bayesian','both')) BayesGLM_left <- BayesGLM_surface(session_data, vertices = verts_left, faces = faces_left, scale=TRUE, num.threads=4, return_INLA_result=FALSE, outfile = NULL) else BayesGLM_left <- NULL

  }


  ### RIGHT HEMISPHERE
  if(do_right){

    #set up mesh
    surf_right <- readGIfTI(fname_gifti_right)$data
    verts_right <- surf_right$pointset
    faces_right <- surf_right$triangle + 1

    #set up session list
    session_data <- vector('list', n_sess)
    names(session_data) <- session_names
    for(ss in 1:n_sess){
      sess <- list(BOLD = t(cifti_right[[ss]]), design=design[[ss]])
      if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
      session_data[[ss]] <- sess
    }

    ### FIT GLM
    if(GLM_method %in% c('classical','both')) classicalGLM_right <- classicalGLM(session_data) else classicalGLM_right <- NULL
    if(GLM_method %in% c('Bayesian','both')) BayesGLM_right <- BayesGLM_surface(session_data, vertices = verts_right, faces = faces_right, scale=TRUE, num.threads=4, return_INLA_result=FALSE, outfile = NULL) else BayesGLM_right <- NULL
  }

  ### SUBCORTICAL
  if(do_sub){

    # Create and Apply Mask

    #include only specified regions
    mask <- array(nifti_labels %in% vol_regions, dim=dim(nifti_labels))
    nvox <- sum(mask)
    ntime <- dim(nifti_data[[1]])[4]
    BOLD_data_vol <- vector('list', n_sess)
    for(ss in 1:n_sess) BOLD_data_vol[[ss]] <- apply(nifti_data[[ss]], 4, function(x, mask){return(x[mask==1])}, mask=mask) #apply mask to each time point (4th dim)

    #identify and remove any zero-variance voxels (important for Bayesian GLM)
    zerovar <- sapply(BOLD_data_vol, rowVars)
    zerovar <- rowSums2(zerovar==0)
    mask[mask==1] <- !zerovar #exclude voxels with zero variance
    nifti_labels <- nifti_labels*mask #this also masks out excluded regions
    nifti_labels_vec <- nifti_labels[mask==TRUE]
    for(ss in 1:n_sess) BOLD_data_vol[[ss]] <- BOLD_data_vol[[ss]][!zerovar,]
    nvox <- sum(mask)

    #set up session list
    session_data <- vector('list', n_sess)
    names(session_data) <- session_names
    for(ss in 1:n_sess){
      sess <- list(BOLD = t(BOLD_data_vol[[ss]]), design=design[[ss]])
      if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
      session_data[[ss]] <- sess
    }

    #set up SPDE
    if(GLM_method %in% c('Bayesian','both')){
      ### SET UP SUBCORTICAL SPDE
    }

    ### FIT GLM
    if(GLM_method %in% c('classical','both')) classicalGLM_vol <- classicalGLM(session_data) else classicalGLM_vol <- NULL

    ### TO DO: Pass through locations, labels & groups_df instead of spde
    if(GLM_method %in% c('Bayesian','both')) BayesGLM_vol <- BayesGLM_vol3D(session_data, spde=spde, scale=TRUE, num.threads=4, return_INLA_result=FALSE, outfile = NULL) else BayesGLM_vol <- NULL
  }

  if(!do_sub) mask <- nifti_labels <- NULL

  classicalGLM_cifti <- BayesGLM_cifti <- vector('list', n_sess)
  names(classicalGLM_cifti) <- names(BayesGLM_cifti) <- session_names
  for(ss in 1:n_sess){
    if(GLM_method %in% c('classical','both')){
      classicalGLM_cifti[[ss]] <- cifti_make(cortex_left = classicalGLM_left$single_session,
                                   cortex_right = classicalGLM_right$single_session,
                                   subcortical = classicalGLM_vol$single_session,
                                   mask = mask,
                                   labels = nifti_labels)
    }
    if(GLM_method %in% c('Bayesian','both')){
      BayesGLM_cifti[[ss]] <- cifti_make(cortex_left = BayesGLM_left$single_session,
                                           cortex_right = BayesGLM_right$single_session,
                                           subcortical = BayesGLM_vol$single_session,
                                           mask = mask,
                                           labels = nifti_labels)
    }
  }

  result <- list(BayesGLM=BayesGLM_cifti, classicalGLM=classicalGLM_cifti)


}



#' Applies spatial Bayesian GLM to task fMRI data on the cortical surface
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @param vertices A Vx3 matrix of vertex locations of the triangular mesh in Euclidean space.
#' @param faces A Wx3 matrix, where each row contains the vertex indices for a given face or triangle in the triangular mesh. W is the number of faces in the mesh.
#' @param mesh A `inla.mesh` object.  Must be provided if and only if `vertices` and `faces` are not.
#' @param mask (Optional) A logical or 0/1 vector of length V indicating which vertices are to be included.
#' @param scale If TRUE, scale timeseries data so estimates represent percent signal change.  Else, do not scale.
#' @param return_INLA_result If TRUE, object returned will include the INLA model object (can be large).  Default is TRUE. Required for running \code{id_activations} on \code{BayesGLM} model object.
#' @param outfile File name where results will be written (for use by \code{BayesGLM_group}).
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
#' @examples \dontrun{}
BayesGLM_surface <- function(data, vertices = NULL, faces = NULL, mesh = NULL, mask = NULL, scale=TRUE, num.threads=4, return_INLA_result=TRUE, outfile = NULL){

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
    warning('No value supplied for outfile, which is required for post-hoc group modeling.')
  }

  if(is.null(mesh)) mesh <- make_mesh(vertices, faces)

  #ID any zero-variance voxels and remove from analysis
  zero_var <- sapply(data, function(x){
    vars <- colVars(x$BOLD)
    return(vars < 1e-6)
    })
  zero_var <- (rowSums(zero_var) > 0) #check whether any vertices have zero variance in any session
  if(sum(zero_var) > 0){
    if(is.null(mask)) mask <- !zero_var else mask <- mask*(!zero_var)
  }

  if(!is.null(mask)) {
    mask <- as.logical(mask)
    mesh <- submesh.mesh(mask, mesh)
    for(s in 1:n_sess){
      data[[s]]$BOLD <- data[[s]]$BOLD[,mask]
    }
    V <- sum(mask)
  }


  # # Mesh masking
  # if(!is.null()){
  #   mesh_orig <- mesh
  #   mesh <- submesh.mesh(mask, mesh)
  #
  #   # # Mask the data too:
  #   # data_orig <- data
  #   # data <-
  # }

  spde <- inla.spde2.matern(mesh)

  #collect data and design matrices
  y_all <- c()
  X_all_list <- NULL

  for(s in 1:n_sess){

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
      names(X_list) <- session_names[s]

      y_all <- c(y_all, y_vec)
      X_all_list <- c(X_all_list, X_list)
  }

  #construct betas and repls objects
  replicates_list <- organize_replicates(n_sess=n_sess, n_task=K, mesh=mesh)
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

  model_data <- make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)

  #estimate model using INLA
  INLA_result <- estimate_model(formula=formula, data=model_data, A=model_data$X, spde=spde, prec_initial=1, num.threads=num.threads)

  #extract useful stuff from INLA model result
  beta_estimates <- extract_estimates(object=INLA_result, session_names=session_names) #posterior means of latent task field
  theta_posteriors <- get_posterior_densities(object=INLA_result, spde) #hyperparameter posterior densities

  #extract stuff needed for group analysis
  mu.theta <- INLA_result$misc$theta.mode
  Q.theta <- solve(INLA_result$misc$cov.intern)


  #construct object to be returned
  if(return_INLA_result){
    result <- list(INLA_result = INLA_result,
                   mesh = mesh,
                   session_names = session_names,
                   beta_names = beta_names,
                   beta_estimates = beta_estimates,
                   theta_posteriors = theta_posteriors,
                   mu.theta = mu.theta, #for joint group model
                   Q.theta = Q.theta, #for joint group model
                   y = y_all, #for joint group model
                   X = X_all_list, #for joint group model
                   call = match.call())
  } else {
    result <- list(INLA_result = NULL,
                   mesh = mesh,
                   session_names = session_names,
                   beta_names = beta_names,
                   beta_estimates = beta_estimates,
                   theta_posteriors = theta_posteriors,
                   mu.theta = mu.theta, #for joint group model
                   Q.theta = Q.theta, #for joint group model
                   y = y_all, #for joint group model
                   X = X_all_list, #for joint group model
                   call = match.call())
  }


  class(result) <- "BayesGLM"

  if(!is.null(outfile)){
    save(outfile, result)
  }

  return(result)

}
