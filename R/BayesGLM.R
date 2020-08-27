#' Performs whole-brain spatial Bayesian GLM for fMRI task activation
#'
#' @param cifti_fname File path (or vector thereof, for multiple-session modeling) of CIFTI-format fMRI timeseries data (*.dtseries.nii).
#' @param surfL_fname File path of GIFTI-format left cortical surface (*.surf.gii). Must be provided if brainstructures includes "left" and GLM_method is "Bayesian" or "both".
#' @param surfR_fname File path of GIFTI-format right cortical surface (*.surf.gii). Must be provided if brainstructures includes "right" and GLM_method is "Bayesian" or "both".
#' @param sphereL_fname File path of GIFTI-format left spherical surface (*.surf.gii) to use for resampling cifti data and gifti surfaces to lower resolution. Must be provided if GLM_method is "Bayesian" or "both" and resample is not NULL.
#' @param sphereR_fname File path of GIFTI-format right spherical surface (*.surf.gii) to use for resampling cifti data and gifti surfaces to lower resolution. Must be provided if GLM_method is "Bayesian" or "both" and resample is not NULL.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param wb_path (Optional) Path to Connectome Workbench folder or executable.
#'  If not provided, should be set with
#'  \code{ciftiTools.setOption("wb_path", "path/to/workbench")}.
#' @param design A TxK task design matrix (or list of such matrices, for multiple-session modeling) with column names representing tasks. Each column represents the expected BOLD response due to each task, a convolution of the hemodynamic response function (HRF) and the task stimulus.  Must be provided if and only if onsets=NULL.  Note that the scale of the regressors will affect the scale and interpretation of the beta coefficients, so imposing a proper scale (e.g., set maximum to 1) is recommended.
#' @param onsets A matrix of onsets (first column) and durations (second column) for each task in SECONDS, organized as a list where each element of the list corresponds to one task. Names of list should be task names. (Or for multi-session modeling, a list of such lists.)  Must be provided if and only if design=NULL.
#' @param TR The temporal resolution of the data in seconds.  Must be provided if onsets provided.
#' @param nuisance (Optional) A TxJ matrix of nuisance signals (or list of such matrices, for multiple-session modeling).
#' @param nuisance_include (Optional) Additional nuisance covariates to include.  Default is 'drift' (linear and quadratic drift terms) and 'dHRF' (temporal derivative of each column of design matrix).
#' @param scale_BOLD If TRUE (default), scale timeseries data so estimates represent percent signal change.  Else, center but do not scale.
#' @param scale_design If TRUE (default), scale design matrix so maximum value is equal to 1.  Else, do not scale.
#' @param GLM_method Either 'Bayesian' for spatial Bayesian GLM only, 'classical' for the classical GLM only, or 'both' to return both classical and Bayesian estimates of task activation.
#' @param session_names (Optional) A vector of names corresponding to each session.
#' @param resamp_res The number of vertices to which each cortical surface should be resampled, or NULL if no resampling is to be performed. For computational feasibility, a value of 10000 or lower is recommended.
#' @param num.threads Maximum number of threads the inla-program will use for model estimation
#' @param verbose Logical indicating if INLA should run in a verbose mode (default FALSE).
#' @param write_dir (Optional) Location where to write resampled data (if resample != NULL) and output files. If NULL, use the current working directory.
#' @param outfile (Optional) File name (without extension) of output file for BayesGLM result to use in Bayesian group modeling.
#' @param return_INLA_result If TRUE, object returned will include the INLA model object (can be large).  Default is FALSE.  Required for running \code{id_activations} on \code{BayesGLM} model object (but not for running \code{BayesGLM_joint} to get posterior quantities of group means or contrasts).
#'
#' @return An object of class BayesGLM, a list containing ...
#' @export
#' @importFrom ciftiTools read_cifti resample_gifti as.xifti
#' @importFrom matrixStats rowVars rowSums2
#' @importFrom gifti readGIfTI
#' @importFrom INLA inla.pardiso.check inla.setOption
#'
#' @details This function uses a system wrapper for the 'wb_command' executable. The user must first download and install the Connectome Workbench,
#' available from https://www.humanconnectome.org/software/get-connectome-workbench. The 'wb_path' argument is the full file path to the 'wb_command' executable file.
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
BayesGLM <- function(cifti_fname,
                     surfL_fname=NULL, surfR_fname=NULL,
                     sphereL_fname=NULL, sphereR_fname=NULL,
                     brainstructures=c('left','right','subcortical'),
                     wb_path=NULL,
                     design=NULL, onsets=NULL, TR=NULL,
                     nuisance=NULL, nuisance_include=c('drift','dHRF'),
                     scale_BOLD=TRUE, scale_design=TRUE,
                     GLM_method='both',
                     session_names=NULL,
                     resamp_res=10000,
                     num.threads=4,
                     verbose=FALSE,
                     write_dir=NULL,
                     outfile=NULL,
                     return_INLA_result=FALSE){

  do_Bayesian <- (GLM_method %in% c('both','Bayesian'))
  do_classical <- (GLM_method %in% c('both','classical'))

  if(do_Bayesian){
    flag <- inla.pardiso.check()
    if(grepl('FAILURE',flag)) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO is required for computational efficiency. See inla.pardiso().')
    inla.setOption(smtp='pardiso')
  }

  if(is.null(write_dir)){
    write_dir <- getwd()
  } else if(!file.exists(write_dir)){
    stop('write_dir does not exist, check and try again.')
  }
  #TO DO: Check that the user has write permissions in write_dir
  #TO DO: Damon: We could integrate ciftiTools::format_path() with this package too.
  # It appends a directory to a file name and checks if it exists/is writeable/is readable.
  #TO DO: Test this line of code for Windows systems.  May need to try both options for winslash argument with a tryCatch or if statement to keep the one that works.
  write_dir <- normalizePath(write_dir) #generate full path

  # Check that arguments are compatible
  brainstructures <- ciftiTools:::match_input(
    brainstructures, c("left","right","subcortical","all"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }
  do_left <- ('left' %in% brainstructures)
  do_right <- ('right' %in% brainstructures)
  do_sub <- ('subcortical' %in% brainstructures)

  if(do_left & is.null(surfL_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  if(do_right & is.null(surfR_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  if((is.null(design) + is.null(onsets)) != 1) stop('design OR onsets must be provided, but not both')
  if(!is.null(onsets) & is.null(TR)) stop('Please provide TR if onsets provided')
  # if(do_sub){
  #   if(length(unique(vol_regions)) != length(vol_regions)) stop('vol_regions must contain no repeated values.')
  #   if(min(is.element(vol_regions, 3:21))==0) stop('vol_regions must include only integer values between 3 and 21.')
  # }

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

  # ### First, perform resampling
  # if(do_Bayesian){
  #   if(!is.null(resamp_res)){
  #     if(!is.numeric(resamp_res)) stop('resamp_res must be numeric')
  #     if(round(resamp_res) != resamp_res) stop('resamp_res must be an integer')
  #     if(resamp_res > 30000 | resamp_res < 1000) stop('resamp_res must be a number between 1,000 and 30,000')
  #     if(is.null(sphereL_fname) | is.null(sphereR_fname)) stop('Must provide sphereL_fname and sphereR_fname to resamp_res.')
  #
  #     #cifti file resampling
  #     cat('\n RESAMPLING CIFTI TIMESERIES FILES TO ', resamp_res, ' RESOLUTION \n')
  #     cifti_fname2 <- rep(NA, n_sess)
  #     cifti_dir <- dirname(cifti_fname[1])
  #     for(ss in 1:n_sess){
  #       cifti_extn <- get_cifti_extn(cifti_fname[ss])
  #       cifti_fname2[ss] <- paste0(gsub(cifti_extn, '', cifti_fname[ss]), resamp_res, '.', cifti_extn)
  #       cifti_fname2[ss] <- file.path(write_dir,basename(cifti_fname2[ss]))
  #       delete_helper_files <- FALSE
  #       if(ss==1){ make_helper_files <- TRUE } else { make_helper_files <- FALSE }
  #       resample_cifti(cifti_orig = cifti_fname[ss],
  #                      cifti_target = cifti_fname2[ss],
  #                      sphere_orig_L = sphereL_fname,
  #                      sphere_orig_R = sphereR_fname,
  #                      target_res = resamp_res,
  #                      wb_path = wb_path,
  #                      make_helper_files = make_helper_files,
  #                      delete_helper_files = delete_helper_files)
  #     }
  #     cifti_fname <- cifti_fname2
  #
  #     #gifti file resampling
  #     cat('\n RESAMPLING GIFTI SURFACE FILES TO ', resamp_res, ' RESOLUTION \n')
  #     fnames_gifti <- c(surfL_fname, surfR_fname, surf2L_fname, surf2R_fname)
  #     fnames_sphere_orig <- c(sphereL_fname, sphereR_fname, sphereL_fname, sphereR_fname)
  #     fnames_sphere_target <- file.path(cifti_dir, 'helper_files_resampling', c('Sphere.target.L.surf.gii',  'Sphere.target.R.surf.gii', 'Sphere.target.L.surf.gii', 'Sphere.target.R.surf.gii'))
  #     notnull <- c(TRUE, TRUE, !is.null(surf2L_fname), !is.null(surf2R_fname))
  #     fnames_sphere_orig <- fnames_sphere_orig[notnull]
  #     fnames_sphere_target <- fnames_sphere_target[notnull]
  #     fnames_gifti_target <- gsub('surf.gii', paste0(resamp_res,'.surf.gii'), fnames_gifti)
  #     for(gg in 1:length(fnames_gifti)){
  #       resample_gifti(gifti_orig = fnames_gifti[gg],
  #                      gifti_target = fnames_gifti_target[gg],
  #                      sphere_orig = fnames_sphere_orig[gg],
  #                      sphere_target = fnames_sphere_target[gg],
  #                      wb_path = wb_path,
  #                      overwrite = FALSE)
  #     }
  #
  #     #redefine gifti file names to resampled files
  #     surfL_fname <- gsub('surf.gii', paste0(resamp_res,'.surf.gii'), surfL_fname)
  #     surfR_fname <- gsub('surf.gii', paste0(resamp_res,'.surf.gii'), surfR_fname)
  #     if(!is.null(surf2L_fname)) surf2L_fname <- gsub('surf.gii', paste0(resamp_res,'.surf.gii'), surf2L_fname)
  #     if(!is.null(surf2R_fname)) surf2R_fname <- gsub('surf.gii', paste0(resamp_res,'.surf.gii'), surf2R_fname)
  #   }
  # }


  cat('\n SETTING UP DATA \n')

  ### For each session, separate the CIFTI data into left/right/sub and read in files
  if(do_left) cifti_left <- vector('list', n_sess)
  if(do_right) cifti_right <- vector('list', n_sess)
  if(do_sub) nifti_data <- nifti_labels <- vector('list', n_sess)

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
        sphereL_fname=sphereL_fname, sphereR_fname=sphereR_fname,
        write_dir=write_dir,
        wb_path=wb_path
      )
      if(do_left) surf_left <- cifti_ss$surf$cortex_left
      if(do_right) surf_right <- cifti_ss$surf$cortex_right
    } else {
      cifti_ss <- read_cifti(
        cifti_fname[ss],
        brainstructures=brainstructures,
        resamp_res=resamp_res,
        sphereL_fname=sphereL_fname, sphereR_fname=sphereR_fname,
        write_dir=write_dir,
        wb_path=wb_path
      )
    }

    if(do_left) {
      cifti_left[[ss]] <- matrix(NA, nrow=length(cifti_ss$meta$cortex$medial_wall_mask$left), ncol=ncol(cifti_ss$data$cortex_left))
      cifti_left[[ss]][cifti_ss$meta$cortex$medial_wall_mask$left,, drop=FALSE] <- cifti_ss$data$cortex_left
    }
    if(do_right) {
      cifti_right[[ss]] <- matrix(NA, nrow=length(cifti_ss$meta$cortex$medial_wall_mask$right), ncol=ncol(cifti_ss$data$cortex_right))
      cifti_right[[ss]][cifti_ss$meta$cortex$medial_wall_mask$right,, drop=FALSE] <- cifti_ss$data$cortex_right
    }
    #if(do_sub) { nifti_data[[ss]] <- cifti_ss$VOL; ntime <- ncol(cifti_ss$VOL) }
    #if(do_sub & ss==1) nifti_labels[[ss]] <- cifti_ss$LABELS

    if(make_design){
      cat(paste0('    Constructing design matrix for session ', ss, '\n'))
      design[[ss]] <- make_HRFs(onsets[[ss]], TR=TR, duration=ntime)
    }

  }
  # #check that labels are the same across all sessions
  # if(do_sub) {
  #   if(n_sess > 1) {
  #     tmp <- sapply(nifti_labels, function(x) {all.equal(x,nifti_labels[[1]])})
  #     if(min(tmp)==0) stop('Subcortical labels must match across all sessions in cifti data. Check compatibility of cifti files.')
  #   }
  #   nifti_labels <- nifti_labels[[1]]
  # }

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
      scale_design <- F
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


  ### LEFT HEMISPHERE
  if(do_left){

    cat('\n ... LEFT CORTEX \n')

    #set up mesh
    #surf_left <- readGIfTI(surfL_fname)$data
    verts_left <- surf_left$vertices #first surface is used for modeling
    faces_left <- surf_left$faces
    #if(min(faces_left)==0) faces_left <- faces_left + 1

    #set up session list
    session_data <- vector('list', n_sess)
    names(session_data) <- session_names
    for(ss in 1:n_sess){
      sess <- list(BOLD = t(cifti_left[[ss]]), design=design[[ss]])
      if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
      session_data[[ss]] <- sess
    }

    ### FIT GLM(s)

    if(!is.null(outfile)) outfile_left <- paste0(outfile, '_left.Rdata') else outfile_left <- NULL

    if(do_classical) classicalGLM_left <- classicalGLM(session_data,
                                                       scale_BOLD=scale_BOLD,
                                                       scale_design = scale_design)
    if(do_Bayesian) BayesGLM_left <- BayesGLM_surface(session_data,
                                                      vertices = verts_left,
                                                      faces = faces_left,
                                                      scale_BOLD=scale_BOLD,
                                                      scale_design = scale_design,
                                                      num.threads=num.threads,
                                                      return_INLA_result=return_INLA_result,
                                                      outfile = file.path(write_dir,outfile_left),
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

    if(!is.null(outfile)) outfile_right <- paste0(outfile, '_right.Rdata') else outfile_right <- NULL

    if(do_classical) classicalGLM_right <- classicalGLM(session_data,
                                                        scale_BOLD=scale_BOLD,
                                                        scale_design = scale_design)
    if(do_Bayesian) BayesGLM_right <- BayesGLM_surface(session_data,
                                                      vertices = verts_right,
                                                      faces = faces_right,
                                                      scale_BOLD=scale_BOLD,
                                                      scale_design = scale_design,
                                                      num.threads=num.threads,
                                                      return_INLA_result=return_INLA_result,
                                                      outfile = file.path(write_dir,outfile_right),
                                                      verbose=verbose)
  }

  # ### SUBCORTICAL
  # if(do_sub){
  #
  #   cat('\n ... SUBCORTICAL \n')
  #
  #   # Create and Apply Mask
  #
  #   #include only specified regions
  #   mask <- array(nifti_labels %in% vol_regions, dim=dim(nifti_labels))
  #   nvox <- sum(mask)
  #   ntime <- dim(nifti_data[[1]])[4]
  #   BOLD_data_vol <- vector('list', n_sess)
  #   for(ss in 1:n_sess) BOLD_data_vol[[ss]] <- apply(nifti_data[[ss]], 4, function(x, mask){return(x[mask==1])}, mask=mask) #apply mask to each time point (4th dim)
  #
  #   #identify and remove any zero-variance voxels (important for Bayesian GLM)
  #   zerovar <- sapply(BOLD_data_vol, rowVars)
  #   zerovar <- rowSums2(zerovar==0)
  #   mask[mask==1] <- !zerovar #exclude voxels with zero variance
  #   nifti_labels <- nifti_labels*mask #this also masks out excluded regions
  #   nifti_labels_vec <- nifti_labels[mask==TRUE]
  #   for(ss in 1:n_sess) BOLD_data_vol[[ss]] <- BOLD_data_vol[[ss]][!zerovar,]
  #   nvox <- sum(mask)
  #
  #   #set up session list
  #   session_data <- vector('list', n_sess)
  #   names(session_data) <- session_names
  #   for(ss in 1:n_sess){
  #     sess <- list(BOLD = t(BOLD_data_vol[[ss]]), design=design[[ss]])
  #     if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
  #     session_data[[ss]] <- sess
  #   }
  #
  #   #set up SPDE
  #   if(do_Bayesian){
  #     ### SET UP SUBCORTICAL SPDE
  #   }
  #
  #   ### FIT GLM
  #   if(do_classical) classicalGLM_vol <- classicalGLM(session_data) else classicalGLM_vol <- NULL
  #
  #   ### TO DO: Pass through locations, labels & groups_df instead of spde
  #   #if(do_Bayesian) BayesGLM_vol <- BayesGLM_vol3D(session_data, spde=spde, scale_BOLD=TRUE, num.threads=num.threads, return_INLA_result=return_INLA_result, outfile = NULL) else BayesGLM_vol <- NULL
  # }
  #
  # if(!do_sub) mask <- nifti_labels <- NULL

  ### CONSTRUCT BETA ESTIMATES AS CIFTI OBJECTS

  cat('\n PUTTING RESULTS IN CIFTI FORMAT \n')

  classicalGLM_cifti <- BayesGLM_cifti <- vector('list', n_sess)
  names(classicalGLM_cifti) <- names(BayesGLM_cifti) <- session_names
  for(ss in 1:n_sess){
    if(do_classical){
      classicalGLM_cifti[[ss]] <- as.xifti(
        cortexL = classicalGLM_left[[ss]],
        cortexR = classicalGLM_right[[ss]]
        #subcortVol = classicalGLM_vol$single_session,
        #mask = mask,
        #subcortLab = nifti_labels
                                             )
    }
    if(do_Bayesian){
      BayesGLM_cifti[[ss]] <- as.xifti(
        cortexL = BayesGLM_left$beta_estimates[[ss]],
        cortexR = BayesGLM_right$beta_estimates[[ss]]
        #subcortVol = BayesGLM_vol$single_session,
        #mask = mask,
        #subcortLab = nifti_labels
                                         )
    }
  }

  result <- list(betas_Bayesian = BayesGLM_cifti,
                 betas_classical = classicalGLM_cifti,
                 GLMs_Bayesian = list(cortexL = BayesGLM_left,
                                     cortexR = BayesGLM_right,
                                     subcortical = BayesGLM_vol),
                 GLMs_classical = list(cortexL = classicalGLM_left,
                                      cortexR = classicalGLM_right,
                                      subcortical = classicalGLM_vol),
                 design = design)

  cat('\n DONE! \n')

  return(result)
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
#' @param scale_BOLD If TRUE, scale timeseries data so estimates represent percent signal change.  Else, center but do not scale.
#' @param num.threads Maximum number of threads the inla-program will use for model estimation
#' @param return_INLA_result If TRUE, object returned will include the INLA model object (can be large).  Default is FALSE.  Required for running \code{id_activations} on \code{BayesGLM} model object (but not for running \code{BayesGLM_joint} to get posterior quantities of group means or contrasts).
#' @param outfile File name where results will be written (for use by \code{BayesGLM_group}).
#' @param verbose Logical indicating if INLA should run in a verbose mode (default FALSE).
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern inla.pardiso.check inla.setOption
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
BayesGLM_surface <- function(data, vertices = NULL, faces = NULL, mesh = NULL, mask = NULL, scale_BOLD=TRUE, scale_design = TRUE, num.threads=4, return_INLA_result=TRUE, outfile = NULL, verbose=FALSE){

  #check whether data is a list OR a session (for single-session analysis)
  #check whether each element of data is a session (use is.session)
  # V = number of data locations
  # T = length of time series for each session (vector)
  # K = number of unique tasks in all sessions

  # Check to see that the INLA package is installed
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("This function requires the INLA package (see www.r-inla.org/download)")


  # Check to see if PARDISO is installed
  if(!exists("inla.pardiso.check", mode = "function")){
    warning("Please update to the latest version of INLA for full functionality and PARDISO compatibility (see www.r-inla.org/download)")
  }else{
    if(inla.pardiso.check() == "FAILURE: PARDISO IS NOT INSTALLED OR NOT WORKING"){
      warning("Consider enabling PARDISO for faster computation (see inla.pardiso())")}
    else {
      inla.setOption(smtp='pardiso')
    }
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
    message('No value supplied for outfile, which is required for post-hoc group modeling.')
  }

  if(is.null(mesh)) mesh <- make_mesh(vertices, faces)

  #ID any zero-variance voxels and remove from analysis
  zero_var <- sapply(data, function(x){
    x$BOLD[is.na(x$BOLD)] <- 0 #to detect medial wall locations coded as NA
    x$BOLD[is.nan(x$BOLD)] <- 0 #to detect medial wall locations coded as NaN
    vars <- colVars(x$BOLD)
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
  cat('\n ...... estimating model with INLA')
  system.time(INLA_result <- estimate_model(formula=formula, data=model_data, A=model_data$X, spde, prec_initial=1, num.threads=num.threads, verbose=verbose))
  cat('\n ...... model estimation completed')

  #extract useful stuff from INLA model result
  beta_estimates <- extract_estimates(object=INLA_result, session_names=session_names, mask=mask) #posterior means of latent task field
  theta_posteriors <- get_posterior_densities(object=INLA_result, spde, beta_names) #hyperparameter posterior densities

  #extract stuff needed for group analysis
  mu.theta <- INLA_result$misc$theta.mode
  Q.theta <- solve(INLA_result$misc$cov.intern)

  #construct object to be returned
  if(return_INLA_result){
    result <- list(INLA_result = INLA_result,
                   mesh = mesh,
                   mask = mask,
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
                   mask = mask,
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
    save(result, file=outfile)
  }

  return(result)

}

#' Estimate INLA model
#'
#' @param formula Formula to put into inla
#' @param data Dataset
#' @param A Large, sparse observation matrix
#' @param spde The spatial model, an object of class inla.spde
#' @param prec_initial Initial precision
#' @param num.threads Number of threads
#' @param int.strategy INLA strategy for numerical integration.  "eb" (empirical Bayes) is recommended for computational efficiency, or "ccd" for greater accuracy
#' @param verbose Logical indicating if should run in a verbose mode (default FALSE).
#'
#' @return Results from INLA
#' @export
#' @importFrom INLA inla
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
estimate_model <- function(formula, data, A, spde, prec_initial, num.threads=4, int.strategy = "eb", verbose=FALSE){

  result <- inla(formula, data=data, control.predictor=list(A=A, compute = TRUE),
                 verbose = verbose, keep = FALSE, num.threads = num.threads,
                 control.inla = list(strategy = "gaussian", int.strategy = int.strategy),
                 control.family=list(hyper=list(prec=list(initial=prec_initial))),
                 control.compute=list(config=TRUE)) #required for excursions
  return(result)
}


#' Make Formula
#'
#' @param beta_names char vector of the names of each bbeta object in the environment
#' @param repl_names char vector of the names of each replicate object in the environment
#' @param hyper_initial Optional vector of initial values for hyperparameters of each latent field OR a list with each element corresponding to one column of the X matrix
#'
#' @return A formula representing the Bayesian GLM to be passed to `inla()`
#'
#' @importFrom stats as.formula
#'
#'
make_formula <- function(beta_names, repl_names, hyper_initial=NULL){

  # Example:
  # beta_names = bbeta1, bbeta2, ...
  # repl_names = repl1, repl2, ...
  # formula: y ~ -1 + f(bbeta1, model = spde, replicate = repl1) + f(bbeta2, model = spde_sh, replicate = repl2)

  # check length of beta_names, repl_names, hyper_initial

  n_beta <- length(beta_names)

  if(!is.null(hyper_initial)){
    #if hyper_list provided is a vector, repeat it n_beta times as a list
    if(!is.list(hyper_initial)){
      hyper_initial <- rep(list(hyper_initial), n_beta)
    }
    hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')
  } else {
    hyper_vec <- NULL
  }

  formula_vec <- paste0('f(',beta_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
  formula_vec <- c('y ~ -1', formula_vec)
  formula_str <- paste(formula_vec, collapse=' + ')
  return(formula_str)
}

#' Make data list to be passed to \code{estimate_model}
#'
#' @param y Vectorized BOLD data (all voxels, sessions, etc.)
#' @param X List (length = number of sessions) of sparse design matrices size TVxVK from each session, each created using `organize_data()`
#' @param betas List (length = number of tasks) of bbeta objects from organize_replicates
#' @param repls List (length = number of tasks) of repl objects from organize_replicates
#'
#' @return List
#'
#' @importFrom Matrix bdiag
#'
#'
make_data_list <- function(y, X, betas, repls){

  # Check length/dimensions of y, X, elements of betas and repls all match
  n_sess <- length(X)
  nx <- length(betas) #check same as length(repls)
  #figure out nvox
  #check dim(X)
  #check length of betas and repls

  numel <- 1 + length(betas) + length(repls) + 1
  model_data <- vector('list', length=numel)
  names(model_data) <- c('y', 'X', names(betas), names(repls))
  model_data$y <- y
  model_data$X <- bdiag(X) #row/col structure: sess1_beta1, sess1_beta2, sess2_beta1, sess2_beta2, ...

  nbeta <- length(betas)
  for(i in 1:nbeta){
    model_data[[2+i]] <- betas[[i]]
    model_data[[2+nbeta+i]] <- repls[[i]]
  }

  return(model_data)
}


#' Extract Estimates of Activation
#'
#' @description Obtains the posterior mean or other summary statistic for each latent field
#'
#' @param object An object of class ‘"inla"’, a result of a call to inla
#' @param session_names Vector of fMRI session names
#' @param mask (Optional) Original mask applied to data before model fitting
#' @param stat A string representing the posterior summary statistic to be returned
#'
#' @return Estimates from inla model
#'
extract_estimates <- function(object, session_names, mask=NULL, stat='mean'){

  if(class(object) != "inla"){
    stop("Object is not of class 'inla'")
  }

  res.beta <- object$summary.random
  nbeta <- length(res.beta)
  beta_names <- names(res.beta)

  n_sess <- length(session_names)
  n_loc <- length(res.beta[[1]]$mean)/n_sess #number of locations for which beta is estimated
  if(!is.null(mask)) {
    V <- length(mask)
    if(sum(mask) != n_loc) warning('Number of nonzeros in mask does not equal the number of data locations in the model')
  } else {
    V <- n_loc
    mask <- rep(1,V)
  }
  betas <- vector('list', n_sess)
  names(betas) <- session_names

  stat_names <- names(res.beta[[1]])
  if(! (stat %in% stat_names) ) stop(paste0('stat must be one of following: ', paste(stat_names, collapse = ', ')))
  stat_ind <- which(stat_names==stat)


  for(v in 1:n_sess){
    inds_v <- (1:n_loc) + (v-1)*n_loc #indices of beta vector corresponding to session v
    betas_v <- matrix(NA, nrow=n_loc, ncol=nbeta)
    colnames(betas_v) <- beta_names

    for(i in 1:nbeta){
      est_iv <- res.beta[[i]][[stat_ind]][inds_v]
      betas_v[,i] <- est_iv
    }
    betas[[v]] <- matrix(NA, nrow=V, ncol=nbeta)
    betas[[v]][mask==1,] <- betas_v
  }
  return(betas)
}


#' Extracts posterior density estimates for hyperparameters
#'
#' @param object An object of class ‘"inla"’, a result of a call to \code{inla()}
#' @param spde The model used for the latent fields in the \code{inla()} call, an object of class ‘"inla.spde"’
#' @param beta_names (Optional) Descriptive names of model regressors (tasks).
#'
#' @return Long-form data frame containing posterior densities for the hyperparameters associated with each latent field
#'
#' @importFrom INLA inla.spde2.result
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
get_posterior_densities <- function(object, spde, beta_names=NULL){

  beta_names_model <- names(object$summary.random)
  numbeta <- length(beta_names_model)
  if(numbeta != length(beta_names)) {
    warning('Length of beta_names invalid.  Setting to NULL.')
    beta_names <- NULL
  }

  for(b in 1:numbeta){
    name_b = beta_names_model[b]
    result.spde.b = inla.spde2.result(object, name_b, spde)
    # Kappa and Tau
    log_kappa.b = as.data.frame(result.spde.b$marginals.log.kappa$kappa.1)
    log_tau.b = as.data.frame(result.spde.b$marginals.log.tau$tau.1)
    names(log_kappa.b) <- names(log_tau.b) <- c('value','density')
    log_kappa.b$param <- 'log_kappa'
    log_tau.b$param <- 'log_tau'
    df.b <- rbind(log_kappa.b, log_tau.b)
    df.b$beta <- name_b
    if(!is.null(beta_names)) df.b$name <- beta_names[b] else df.b$name <- NA
    if(b == 1) df <- df.b else df <- rbind(df, df.b)
  }
  df <- df[,c('beta','name','param','value','density')]
  return(df)
}


#' Extracts posterior density estimates for hyperparameters
#'
#' @param object An object of class ‘"inla"’, a result of a call to \code{inla()}
#' @param spde The model used for the latent fields in the \code{inla()} call, an object of class ‘"inla.spde"’
#'
#' @return Long-form data frame containing posterior densities for the hyperparameters associated with each latent field
#' @export
#' @importFrom INLA inla.spde2.result
#' @importFrom INLA inla.extract.el
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
get_posterior_densities_vol3D <- function(object, spde){

  hyper_names <- names(object$marginals.hyperpar)

  for(h in hyper_names){
    df.h <- inla.extract.el(object$marginals.hyperpar, h)
    df.h <- as.data.frame(df.h)
    names(df.h) <- c('x','y')
    df.h$hyper_name <- h
    df.h$beta_name <- ifelse(grepl('bbeta',h), #if bbeta appears in name
                             gsub('.+bbeta','bbeta',h), #rename as bbeta*
                             NA) #NA if this is the precision hyperparameter
    df.h$theta_name <- ifelse(grepl('Theta',h), #if bbeta appears in name
                              gsub(' for.+','',h), #rename as Theta*
                              NA) #NA if this is the precision hyperparameter

    if(h==hyper_names[1]) df <- df.h else df <- rbind(df, df.h)
  }
  return(df)
}

