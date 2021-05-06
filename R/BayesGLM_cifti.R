#' BayesGLM for CIFTI
#'
#' Performs spatial Bayesian GLM on the cortical surface for fMRI task activation
#'
#' @section Connectome Workbench Requirement:
#'  This function uses a system wrapper for the 'wb_command' executable. The
#'  user must first download and install the Connectome Workbench, available
#'  from https://www.humanconnectome.org/software/get-connectome-workbench .
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
#' @param design,onsets,TR Either provide \code{design}, or provide both \code{onsets} and \code{TR}.
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
#' @param nuisance_include (Optional) Additional nuisance covariates to include.
#'  Default is 'drift' (linear and quadratic drift terms) and 'dHRF' (temporal
#'  derivative of each column of design matrix). Set to \code{NULL} to not do any
#'  nuisance regressors, or just one of these to use only one.
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @param GLM_method Either 'Bayesian' for spatial Bayesian GLM only, 'classical' for the classical GLM only, or 'both' to return both classical and Bayesian estimates of task activation.
#' @param ar_order (numeric) Controls prewhitening. If greater than zero, this
#'  should be a number indicating the order of the autoregressive model to use
#'  for prewhitening. If zero, do not prewhiten. Default: \code{6}.
#' @param ar_smooth FWHM parameter for smoothing. Remember that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} or \code{NULL}
#'  to not do any smoothing. Default: \code{5}.
#' @param session_names (Optional) A vector of names corresponding to each
#'   session.
#' @param resamp_res The number of vertices to which each cortical surface should be resampled, or NULL if no resampling is to be performed. For computational feasibility, a value of 10000 or lower is recommended.
#' @inheritParams num.threads_Param
#' @inheritParams verbose_Param_inla
#' @param outfile (Optional) File name (without extension) of output file for
#'  \code{"BayesGLM"} result to use in Bayesian group modeling.
#'  \code{"_left.rds"} or \code{"_right.rds"} will be appended for the left
#'  cortex and right cortex results, respectively. Default: \code{NULL}
#'  (do not save the results to any file).
#' @inheritParams return_INLA_result_Param_FALSE
#' @inheritParams avg_sessions_Param
#' @param trim_INLA (logical) should the \code{INLA_result} objects within the
#'   result be trimmed to only what is necessary to use `id_activations()`? Default: `TRUE`.
#'
#' @return An object of class \code{"BayesGLM"}, a list containing...
#'
# @importFrom ciftiTools read_cifti resample_gifti as.xifti remove_xifti
#' @import ciftiTools
#' @importFrom matrixStats rowVars rowSums2 colVars
#' @importFrom INLA inla.pardiso.check inla.setOption
#' @importFrom parallel detectCores
#'
#' @export
BayesGLM_cifti <- function(cifti_fname,
                           surfL_fname=NULL, surfR_fname=NULL,
                           brainstructures=c('left','right'),
                           design=NULL, onsets=NULL, TR=NULL,
                           nuisance=NULL, nuisance_include=c('drift','dHRF'),
                           scale_BOLD=TRUE, scale_design=TRUE,
                           GLM_method='both',
                           ar_order = 6,
                           ar_smooth = 5,
                           session_names=NULL,
                           resamp_res=10000,
                           num.threads=4,
                           verbose=FALSE,
                           outfile=NULL,
                           return_INLA_result=FALSE,
                           avg_sessions = TRUE,
                           trim_INLA = TRUE){

  # allows for list input
  cifti_fname <- as.character(cifti_fname)

  GLM_method = match.arg(GLM_method, c('both','Bayesian','classical'))

  do_Bayesian <- (GLM_method %in% c('both','Bayesian'))
  do_classical <- (GLM_method %in% c('both','classical'))

  check_BayesGLM(require_PARDISO=do_Bayesian)
  prewhiten <- ar_order > 0

  avail_cores <- parallel::detectCores()
  if(avail_cores < 2) {
    num.threads <- 1
  }

  # Check that arguments are compatible
  brainstructures <- ciftiTools:::match_input(
    brainstructures, c("left","right"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }
  #do_left <- ('left' %in% brainstructures)
  #do_right <- ('right' %in% brainstructures)
  #do_sub <- FALSE
  #do_sub <- ('subcortical' %in% brainstructures)
  if(!is.null(onsets)){
    #for multiple session data, onsets is a list (representing sessions) of lists (representing tasks)
    if(class(onsets[[1]]) == 'list') {
      if(is.null(names(onsets[[1]])))
        beta_names <- paste0("beta",seq_len(length(onsets[[1]])))
      if(!is.null(names(onsets[[1]])))
        beta_names <- names(onsets[[1]])
    }
    #for single session data, onsets is a list (representing tasks) of data frames or matrices
    if(('data.frame' %in% class(onsets[[1]])) | ('matrix' %in% class(onsets[[1]]))) {
      if(is.null(names(onsets))) beta_names <- paste0("beta",seq_len(length(onsets)))
      if(!is.null(names(onsets))) beta_names <- names(onsets)
    }
  }
  if(!is.null(design)) {
    #for multi-session data
    if(class(design) == "list") {
      if(is.null(colnames(design[[1]])))
        beta_names <- paste0("beta",seq_len(ncol(design[[1]])))
      if(!is.null(colnames(design[[1]])))
        beta_names <- colnames(design[[1]])
    }
    #for single session data
    if("matrix" %in% class(design) | "data.frame" %in% class(design)) {
      if(is.null(colnames(design))) beta_names <- paste0("beta",seq_len(ncol(design)))
      if(!is.null(colnames(design))) beta_names <- colnames(design)
    }
  }

  if('left' %in% brainstructures & is.null(surfL_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  if('right' %in% brainstructures & is.null(surfR_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  if((is.null(design) + is.null(onsets)) != 1) stop('design OR onsets must be provided, but not both')
  if(!is.null(onsets) & is.null(TR)) stop('Please provide TR if onsets provided')
  # if(do_sub){
  #   if(length(unique(vol_regions)) != length(vol_regions)) stop('vol_regions must contain no repeated values.')
  #   if(min(is.element(vol_regions, 3:21))==0) stop('vol_regions must include only integer values between 3 and 21.')
  # }

  # Name sessions and check compatibility of multi-session arguments
  n_sess <- length(cifti_fname)
  if(n_sess == 1 & avg_sessions) avg_sessions <- FALSE
  if(n_sess==1){
    if(is.null(session_names)) session_names <- 'single_session'
    if(!is.null(design)) design <- list(single_session = design)
    if(!is.null(onsets)) onsets <- list(single_session = onsets)
    if(!is.null(nuisance)) nuisance <- list(single_session = nuisance)
  } else {
    if(is.null(session_names)) session_names <- paste0('session', 1:n_sess)
    if(!is.null(design)){ if(length(design) != n_sess) stop('If multiple sessions provided (because cifti_fname is a vector), design must be a list of length equal to the number of sessions (or NULL, if onsets provided).') }
    if(!is.null(onsets)){ if(length(onsets) != n_sess) stop('If multiple sessions provided (because cifti_fname is a vector), onsets must be a list of length equal to the number of sessions (or NULL, if design provided).') }
    if(!is.null(nuisance)){ if(length(nuisance) != n_sess) stop('If multiple sessions provided (because cifti_fname is a vector), nuisance must be a list of length equal to the number of sessions (or NULL).') }
  }
  if(length(session_names) != n_sess) stop('If session_names is provided, it must be of the same length as cifti_fname')

  cat('\n SETTING UP DATA \n')

  ### For each session, separate the CIFTI data into left/right/sub and read in files
  #if(do_left) cifti_left <- vector('list', n_sess)
  #if(do_right) cifti_right <- vector('list', n_sess)
  #if(do_sub) nifti_data <- nifti_labels <- vector('list', n_sess)

  # if(is.null(design)) {
  #   make_design <- TRUE
  #   design <- vector('list', length=n_sess)
  # }

  # for(ss in 1:n_sess){

  #   cat(paste0('    Reading in data for session ', ss,'\n'))

  #   if(ss==1){
  #     cifti_ss <- read_cifti(
  #       cifti_fname[ss],
  #       surfL_fname=surfL_fname, surfR_fname=surfR_fname,
  #       brainstructures=brainstructures,
  #       resamp_res=resamp_res
  #     )
  #     if(do_left) surf_left <- cifti_ss$surf$cortex_left
  #     if(do_right) surf_right <- cifti_ss$surf$cortex_right
  #   } else {
  #     cifti_ss <- read_cifti(
  #       cifti_fname[ss],
  #       brainstructures=brainstructures,
  #       resamp_res=resamp_res
  #     )
  #   }

  #   if(do_left) {
  #     cifti_left[[ss]] <- matrix(NA, nrow=length(cifti_ss$meta$cortex$medial_wall_mask$left), ncol=ncol(cifti_ss$data$cortex_left))
  #     cifti_left[[ss]][cifti_ss$meta$cortex$medial_wall_mask$left,] <- cifti_ss$data$cortex_left
  #     ntime <- ncol(cifti_left[[ss]])
  #   }
  #   if(do_right) {
  #     cifti_right[[ss]] <- matrix(NA, nrow=length(cifti_ss$meta$cortex$medial_wall_mask$right), ncol=ncol(cifti_ss$data$cortex_right))
  #     cifti_right[[ss]][cifti_ss$meta$cortex$medial_wall_mask$right,] <- cifti_ss$data$cortex_right
      # ntime <- ncol(cifti_right[[ss]])
    # }
    #if(do_sub) { nifti_data[[ss]] <- cifti_ss$VOL; ntime <- ncol(cifti_ss$VOL) }
    #if(do_sub & ss==1) nifti_labels[[ss]] <- cifti_ss$LABELS

#    if(make_design){
#   cat(paste0('    Constructing design matrix for session ', ss, '\n'))
#    design[[ss]] <- make_HRFs(onsets[[ss]], TR=TR, duration=ntime)
#   }

cat("READING IN CIFTI DATA \n")
cifti_data <- sapply(brainstructures, function(each_hem) {
  sapply(1:n_sess, function(ss){
    cifti_ss <- read_cifti(
      cifti_fname[ss],
      surfL_fname=surfL_fname, surfR_fname=surfR_fname,
      brainstructures=brainstructures,
      resamp_res=resamp_res
    )
    mwall_mask <- cifti_ss$meta$cortex$medial_wall_mask[[each_hem]]
    cdata <- cifti_ss$data[[paste0("cortex_",each_hem)]]
    cifti_out <- matrix(NA,
                        nrow = length(mwall_mask),
                        ncol=ncol(cdata))
    cifti_out[mwall_mask,] <- cdata
    return(list(
      data = cifti_out,
      surf = cifti_ss$surf[[paste0("cortex_", each_hem)]],
      ntime = ncol(cdata),
      cifti = cifti_ss
    ))
  }, simplify = F)
}, simplify = F)

ntimes <- sapply(cifti_data, function(h) sapply(h, `[[`, i = 3), simplify = F)

cat("MAKING DESIGN MATRICES \n")
design <- mapply(make_HRFs, onsets, TR = TR, duration = ntimes[[1]], SIMPLIFY = F)

  #}
  # #check that labels are the same across all sessions (subcortical)
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
    if(length(beta_names) == 1) {
      num_names <- length(unique(tmp))
      if(num_names > 1) stop('task names must match across sessions for multi-session modeling')
    } else {
      num_names <- apply(tmp, 1, function(x) length(unique(x))) #number of unique names per row
      if(max(num_names) > 1) stop('task names must match across sessions for multi-session modeling')
    }
  }

  cat('\n RUNNING MODELS \n')
  classicalGLM_results <- list(left = NULL, right = NULL, vol = NULL)
  BayesGLM_results <- list(left = NULL, right = NULL, vol = NULL)

  # classicalGLM_left <- classicalGLM_right <- classicalGLM_vol <- NULL
  # BayesGLM_left <- BayesGLM_right <- BayesGLM_vol <- NULL

  ### FORMAT DESIGN MATRIX
  # for(ss in 1:n_sess){
  #   if(scale_design){
  #     design[[ss]] <- scale_design_mat(design[[ss]])
  #   } else {
  #     design[[ss]] <- scale(design[[ss]], scale=FALSE) #center design matrix to eliminate baseline
  #   }
  # }
  if(scale_design) design <- sapply(design,scale_design_mat, simplify = F)
  if(!scale_design) design <- sapply(design,scale, scale = F, simplify = F)

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

  #scale_design <- FALSE # This is done to prevent double-scaling in the BayesGLM function

  # >> Loop through brainstructures to complete the analyses on the different hemispheres ----
for(each_hem in brainstructures) {
  cat("\n ...",toupper(each_hem),"CORTEX ANALYSIS \n")
  verts <- cifti_data[[each_hem]][[1]]$surf$vertices
  faces <- cifti_data[[each_hem]][[1]]$surf$faces

  #set up session list
  session_data <- vector('list', n_sess)
  names(session_data) <- session_names
  for(ss in 1:n_sess){
    sess <- list(BOLD = t(cifti_data[[each_hem]][[ss]]$data), design=design[[ss]])
    if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
    session_data[[ss]] <- sess
  }

  if(!is.null(outfile)) {
    if (endsWith(outfile, ".rds")) {
      outfile_name <- gsub(".rds$", paste0("_",each_hem,".rds"), outfile)
    } else {
      outfile_name <- paste0(outfile, "_",each_hem,".rds")
    }
  } else {
    outfile_name <- NULL
  }

  scale_BOLD_hem <- scale_BOLD
  other_hem <- grep(each_hem, c("left","right"), value = T, invert = T)
  # if(prewhiten) {
  #   pw_data_hem <- prewhiten_cifti(data = session_data,
  #                                   scale_BOLD = scale_BOLD_hem,
  #                                   scale_design = FALSE,
  #                                   ar_order = ar_order,
  #                                   ar_smooth =  ar_smooth,
  #                                   hemisphere = each_hem,
  #                                   cifti_data = ciftiTools::remove_xifti(cifti_data[[each_hem]][[1]]$cifti, paste0("cortex_",other_hem)),
  #                                   num.threads = num.threads)
  #   scale_BOLD_hem <- FALSE # done in prewhitening
  #   session_data <- pw_data_hem$data
  # }
  # if(do_classical) classicalGLM_out <- classicalGLM(data=session_data,
  #                                                    scale_BOLD=scale_BOLD_hem,
  #                                                    scale_design = FALSE) # done above
  mask_hem <- make_mask(session_data)
  BayesGLM_out <- BayesGLM(data = session_data,
                          beta_names = beta_names,
                          vertices = verts,
                          faces = faces,
                          mask = mask_hem,
                          scale_BOLD = scale_BOLD_hem,
                          scale_design = FALSE, # done above
                          GLM_method = GLM_method,
                          ar_order = ar_order,
                          ar_smooth = ar_smooth,
                          num.threads = num.threads,
                          return_INLA_result = return_INLA_result,
                          outfile = outfile_name,
                          verbose = verbose,
                          avg_sessions = avg_sessions,
                          trim_INLA = trim_INLA)
  BayesGLM_results[[each_hem]] <- BayesGLM_out[-grep("classical", names(BayesGLM_out))]
  classicalGLM_results[[each_hem]] <- BayesGLM_out[[grep("classical", names(BayesGLM_out))]]
}

  # ### >> LEFT HEMISPHERE ----
  # if(do_left){

  #   cat('\n ... LEFT CORTEX \n')

  #   #set up mesh
  #   #surf_left <- readGIfTI(surfL_fname)$data
  #   verts_left <- surf_left$vertices #first surface is used for modeling
  #   faces_left <- surf_left$faces
  #   #if(min(faces_left)==0) faces_left <- faces_left + 1

  #   #set up session list
  #   session_data <- vector('list', n_sess)
  #   names(session_data) <- session_names
  #   for(ss in 1:n_sess){
  #     sess <- list(BOLD = t(cifti_left[[ss]]), design=design[[ss]])
  #     if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
  #     session_data[[ss]] <- sess
  #   }

  #   ### CHECK FOR FLAT/NA/NaN VERTICES

  #   mask_left <- make_mask(data = session_data)

  #   ### FIT GLM(s)

  #   if(!is.null(outfile)) {
  #     if (endsWith(outfile, ".rds")) {
  #       outfile_left <- gsub(".rds$", "_left.rds", outfile)
  #     } else {
  #       outfile_left <- paste0(outfile, "_left.rds")
  #     }
  #   } else {
  #     outfile_left <- NULL
  #   }

  #   # scale_BOLD_left <- scale_BOLD
  #   # if(prewhiten) {
  #   #   pw_data_left <- prewhiten_cifti(data = session_data,
  #   #                                   scale_BOLD = scale_BOLD_left,
  #   #                                   scale_design = FALSE,
  #   #                                   ar_order = ar_order,
  #   #                                   ar_smooth =  ar_smooth,
  #   #                                   hemisphere = 'left',
  #   #                                   cifti_data = ciftiTools::remove_xifti(cifti_ss, "cortex_right"),
  #   #                                   num.threads = num.threads)
  #   #   scale_BOLD_left <- FALSE # done in prewhitening
  #   #   session_data <- pw_data_left$data
  #   # }


  #   if(do_Bayesian) {
  #     #BayesGLM performs Bayesian and classical GLM
  #     BayesGLM_left <- BayesGLM(data = session_data,
  #                               beta_names = beta_names,
  #                               vertices = verts_left,
  #                               faces = faces_left,
  #                               mask = mask_left,
  #                               scale_BOLD = scale_BOLD,
  #                               scale_design = FALSE, # done above
  #                               ar_order = ar_order,
  #                               ar_smooth = ar_smooth,
  #                               num.threads = num.threads,
  #                               return_INLA_result = return_INLA_result,
  #                               outfile = outfile_left,
  #                               verbose = verbose,
  #                               avg_sessions = avg_sessions,
  #                               trim_INLA = trim_INLA)
  #     classicalGLM_left <- BayesGLM_left$result_classical
  #   }

  #   if(GLM_method=='classical') classicalGLM_left <- classicalGLM(data=session_data,
  #                                                                 beta_names = beta_names,
  #                                                                 vertices = verts_left,
  #                                                                 faces = faces_left,
  #                                                                 mask = mask_left,
  #                                                                 scale_BOLD=scale_BOLD,
  #                                                                 # These will need to be uncommented once the classicalGLM function is complete
  #                                                                 # ar_order=ar_order,
  #                                                                 # ar_smooth=ar_smooth,
  #                                                                 scale_design = FALSE) # done above

  # }

  # ### >> RIGHT HEMISPHERE ----
  # if(do_right){

  #   cat('\n ... RIGHT CORTEX \n')

  #   #set up mesh
  #   verts_right <- surf_right$vertices #first surface is used for modeling
  #   faces_right <- surf_right$faces

  #   #set up session list
  #   session_data <- vector('list', n_sess)
  #   names(session_data) <- session_names
  #   for(ss in 1:n_sess){
  #     sess <- list(BOLD = t(cifti_right[[ss]]), design=design[[ss]])
  #     if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
  #     session_data[[ss]] <- sess
  #   }

  #   ### CHECK FOR FLAT/NA/NaN VERTICES

  #   mask_right <- make_mask(data = session_data)

  #   ### FIT GLM(s)

  #   if(!is.null(outfile)) {
  #     if (endsWith(outfile, ".rds")) {
  #       outfile_right <- gsub(".rds$", "_right.rds", outfile)
  #     } else {
  #       outfile_right <- paste0(outfile, "_right.rds")
  #     }
  #   } else {
  #     outfile_right <- NULL
  #   }

  #   if(do_Bayesian) {
  #     BayesGLM_right <- BayesGLM(data = session_data,
  #                                beta_names = beta_names,
  #                                vertices = verts_right,
  #                                faces = faces_right,
  #                                mask = mask_right,
  #                                scale_BOLD = scale_BOLD,
  #                                scale_design = FALSE, # done above
  #                                ar_order = ar_order,
  #                                ar_smooth = ar_smooth,
  #                                num.threads = num.threads,
  #                                return_INLA_result = return_INLA_result,
  #                                outfile = outfile_right,
  #                                verbose = verbose,
  #                                avg_sessions = avg_sessions,
  #                                trim_INLA = trim_INLA)
  #     classicalGLM_right <- BayesGLM_right$result_classical
  #   }

  #   if(GLM_method=='classical') classicalGLM_right <- classicalGLM(data=session_data,
  #                                                                  scale_BOLD=scale_BOLD,
  #                                                                  # These will need to be uncommented once the classicalGLM function is complete
  #                                                                  # ar_order=ar_order,
  #                                                                  # ar_smooth=ar_smooth,
  #                                                                  scale_design = FALSE) # done above

  # }

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
  cifti_ss <- sapply(cifti_data, function(x) return(x[[1]]$cifti), simplify = F)
  if(length(cifti_ss) == 1) cifti_ss <- cifti_ss[[1]]
  classicalGLM_cifti <- BayesGLM_cifti <- vector('list', n_sess)
  names(classicalGLM_cifti) <- names(BayesGLM_cifti) <- session_names
  datL <- datR <- cortexL_mwall <- cortexR_mwall <- NULL
  if('left' %in% brainstructures) cortexL_mwall <- as.numeric(cifti_ss$meta$cortex$medial_wall_mask$left)
  if('right' %in% brainstructures) cortexR_mwall <- as.numeric(cifti_ss$meta$cortex$medial_wall_mask$right)
  for(ss in 1:n_sess){
    if(do_classical){
      if('left' %in% brainstructures) datL <- classicalGLM_results$left[[ss]]$estimates[cortexL_mwall==1,]
      if('right' %in% brainstructures) datR <- classicalGLM_results$right[[ss]]$estimates[cortexR_mwall==1,]
      classicalGLM_cifti[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = cortexL_mwall,
        cortexR = datR,
        cortexR_mwall = cortexR_mwall
        #subcortVol = classicalGLM_vol$single_session,
        #mask = mask,
        #subcortLab = nifti_labels
      )

      classicalGLM_cifti[[ss]]$meta$cifti$names <- beta_names
    }
    if(do_Bayesian){
      if('left' %in% brainstructures) datL <- BayesGLM_results$left$beta_estimates[[ss]][cortexL_mwall==1,]
      if('right' %in% brainstructures) datR <- BayesGLM_results$right$beta_estimates[[ss]][cortexR_mwall==1,]
      BayesGLM_cifti[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = cortexL_mwall,
        cortexR = datR,
        cortexR_mwall = cortexR_mwall
        #subcortVol = BayesGLM_vol$single_session,
        #mask = mask,
        #subcortLab = nifti_labels
      )
      BayesGLM_cifti[[ss]]$meta$cifti$names <- beta_names
    }
  }

  if(avg_sessions) {
    if(do_classical) {
      # Need to include the missing locations for medial walls within the
      # linear combinations or the xifti object will be mis-mapped.
      if('left' %in% brainstructures) {
        datL <- classicalGLM_results$left$avg$estimates[cortexL_mwall==1,]
        # datL <- matrix(NA, sum(cortexL_mwall),length(beta_names))
        # maskL <- classicalGLM_left$avg$mask[cortexL_mwall == 1]
        # datL[maskL,] <- classicalGLM_left$avg$estimates[cortexL_mwall == 1,]
      }
      if('right' %in% brainstructures) {
        datR <- classicalGLM_results$right$avg$estimates[cortexR_mwall==1,]
        # datR <- matrix(NA, sum(cortexR_mwall),length(beta_names))
        # maskR <- classicalGLM_right$avg$mask[cortexR_mwall == 1]
        # datR[maskR,] <- classicalGLM_right$avg$estimates[cortexR_mwall == 1,]
      }
      # Adding the averages to the front of the BayesGLM_cifti object
      classicalGLM_cifti <- c(
        list(avg = as.xifti(
          cortexL = datL,
          cortexL_mwall = cortexL_mwall,
          cortexR = datR,
          cortexR_mwall = cortexR_mwall
        )),
        classicalGLM_cifti
      )
      classicalGLM_cifti[[1]]$meta$cifti$names <- beta_names
    }

    if (do_Bayesian) {
      # Need to include the missing locations for medial walls within the
      # linear combinations or the xifti object will be mis-mapped.
      if('left' %in% brainstructures) {
        datL <- matrix(NA, sum(cortexL_mwall),length(beta_names))
        maskL <- BayesGLM_results$left$mask[cortexL_mwall == 1]
        datL[maskL,] <- BayesGLM_results$left$avg_beta_estimates
      }
      if('right' %in% brainstructures) {
        datR <- matrix(NA, sum(cortexR_mwall),length(beta_names))
        maskR <- BayesGLM_results$right$mask[cortexR_mwall == 1]
        datR[maskR,] <- BayesGLM_results$right$avg_beta_estimates
      }
      # Adding the averages to the front of the BayesGLM_cifti object
      BayesGLM_cifti <- c(
        list(avg = as.xifti(
          cortexL = datL,
          cortexL_mwall = cortexL_mwall,
          cortexR = datR,
          cortexR_mwall = cortexR_mwall
        )),
        BayesGLM_cifti
      )
      BayesGLM_cifti[[1]]$meta$cifti$names <- beta_names
    }
  }


  # prewhitening_info <- list()
  # if(prewhiten) {
  #   # I take out the first element here because the first argument
  #   # is a copy of the data, which are already provided.
  #   if(do_left) prewhitening_info$left <- pw_data_left[-1]
  #   if(do_right) prewhitening_info$right <- pw_data_right[-1]
  # }

  result <- list(session_names = session_names,
                 beta_names = beta_names,
                 betas_Bayesian = BayesGLM_cifti,
                 betas_classical = classicalGLM_cifti,
                 GLMs_Bayesian = list(cortexL = BayesGLM_results$left,
                                      cortexR = BayesGLM_results$right),
                 #subcortical = BayesGLM_vol),
                 GLMs_classical = list(cortexL = classicalGLM_results$left,
                                       cortexR = classicalGLM_results$right),
                 #subcortical = classicalGLM_vol),
                 #prewhitening_info = prewhitening_info,
                 design = design)

  cat('\n DONE! \n')

  class(result) <- "BayesGLM_cifti"
  return(result)
}
