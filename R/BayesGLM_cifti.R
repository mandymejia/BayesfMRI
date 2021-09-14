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
#' @param cifti_fname File path (or vector thereof, for multiple sessions) of CIFTI-format fMRI timeseries data (*.dtseries.nii).
#' @param surfL_fname File path of GIFTI-format left cortical surface (*.surf.gii). Must be provided if brainstructures includes "left" and Bayes=TRUE.
#' @param surfR_fname File path of GIFTI-format right cortical surface (*.surf.gii). Must be provided if brainstructures includes "right" and Bayes=TRUE.
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
#' @inheritParams Bayes_Param
#' @param ar_order (numeric) Controls prewhitening. If greater than zero, this
#'  should be a number indicating the order of the autoregressive model to use
#'  for prewhitening. If zero, do not prewhiten. Default: \code{6}.
#' @param ar_smooth FWHM parameter for smoothing. Remember that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} or \code{NULL}
#'  to not do any smoothing. Default: \code{5}.
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
#' @param session_names (Optional) A vector of names corresponding to each
#'   session. Ignored if \code{avg_sessions == TRUE}.
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
                           surfL_fname=NULL,
                           surfR_fname=NULL,
                           brainstructures=c('left','right'),
                           design=NULL,
                           onsets=NULL,
                           TR=NULL,
                           nuisance=NULL,
                           nuisance_include=c('drift','dHRF'),
                           scale_BOLD=TRUE,
                           scale_design=TRUE,
                           Bayes=TRUE,
                           ar_order = 6,
                           ar_smooth = 5,
                           resamp_res=10000,
                           num.threads=4,
                           verbose=FALSE,
                           outfile=NULL,
                           return_INLA_result=FALSE,
                           avg_sessions = TRUE,
                           session_names=NULL,
                           trim_INLA = TRUE){

  # allows for list input
  cifti_fname <- as.character(cifti_fname)

  do_Bayesian <- as.logical(Bayes)

  check_BayesGLM(require_PARDISO=do_Bayesian)
  prewhiten <- ar_order > 0

  avail_cores <- parallel::detectCores()
  num.threads <- min(num.threads, avail_cores)
  #if(avail_cores < 2) num.threads <- 1

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
    if('list' %in% class(design)) {
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

  if(do_left & is.null(surfL_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  if(do_right & is.null(surfR_fname)) stop('surfR_fname must be provided if brainstructures includes "right"')

  if((is.null(design) + is.null(onsets)) != 1) stop('design OR onsets must be provided, but not both')
  if(!is.null(onsets) & is.null(TR)) stop('Please provide TR if onsets provided')

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
  if(do_left) BOLD_left <- vector('list', n_sess) else BOLD_left <- NULL
  if(do_right) BOLD_right <- vector('list', n_sess) else BOLD_right <- NULL

  for(ss in 1:n_sess){

    if(n_sess > 1) cat(paste0(' .. reading in data for session ', ss,'\n'))

    if(ss==1){
      #only read in surfaces with first session
      cifti_ss <- read_cifti(
        cifti_fname[ss],
        surfL_fname=surfL_fname, surfR_fname=surfR_fname,
        brainstructures=brainstructures,
        resamp_res=resamp_res
      )
      #grab surface geometry data
      if(do_left) surf_left <- cifti_ss$surf$cortex_left else surf_left <- NULL
      if(do_right) surf_right <- cifti_ss$surf$cortex_right else surf_right <- NULL
      #grab medial wall masks
      if(do_left) cortexL_mwall <- as.numeric(cifti_ss$meta$cortex$medial_wall_mask$left) else cortexL_mwall <- NULL
      if(do_right) cortexR_mwall <- as.numeric(cifti_ss$meta$cortex$medial_wall_mask$right) else cortexR_mwall <- NULL
    } else {
      cifti_ss <- read_cifti(
        cifti_fname[ss],
        brainstructures=brainstructures,
        resamp_res=resamp_res
      )
      #TO DO: check here that medial wall matches first session?
    }

    #grab BOLD data (input NAs in medial wall locations)
    if(do_left) {
      BOLD_left[[ss]] <- matrix(NA, nrow=length(cifti_ss$meta$cortex$medial_wall_mask$left), ncol=ncol(cifti_ss$data$cortex_left))
      BOLD_left[[ss]][cifti_ss$meta$cortex$medial_wall_mask$left,] <- cifti_ss$data$cortex_left
      ntime <- ncol(BOLD_left[[ss]])
    }
    if(do_right) {
      BOLD_right[[ss]] <- matrix(NA, nrow=length(cifti_ss$meta$cortex$medial_wall_mask$right), ncol=ncol(cifti_ss$data$cortex_right))
      BOLD_right[[ss]][cifti_ss$meta$cortex$medial_wall_mask$right,] <- cifti_ss$data$cortex_right
      ntime <- ncol(BOLD_right[[ss]])
    }
  }

  #put together for looping over hemispheres
  BOLD_list <- surf_list <- list(left = NULL, right = NULL)
  BOLD_list$left <- BOLD_left #if RHS if NULL, 'left' will be removed from list
  BOLD_list$right <- BOLD_right #if RHS if NULL, 'right' will be removed from list
  surf_list$left <- surf_left #if RHS if NULL, 'left' will be removed from list
  surf_list$right <- surf_right #if RHS if NULL, 'right' will be removed from list
  rm(BOLD_left, BOLD_right, surf_left, surf_right); gc()

  if(is.null(design)) {
    cat(" MAKING DESIGN MATRICES \n")
    #currently assumes that all sessions have the same duration, could relax this
    #ntimes <- sapply(cifti_data, function(h) sapply(h, `[[`, i = 3), simplify = F)
    design <- mapply(make_HRFs, onsets, TR = TR, duration = ntime, SIMPLIFY = F)
  }

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

  cat(' RUNNING MODELS \n')
  classicalGLM_results <- list(left = NULL, right = NULL)
  BayesGLM_results <- list(left = NULL, right = NULL)

  if(scale_design) design <- sapply(design, scale_design_mat, simplify = F)
  if(!scale_design) design <- sapply(design, scale, scale = F, simplify = F) #center but do not scale

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

  # >> Loop through brainstructures to complete the analyses on the different hemispheres ----
  for(each_hem in brainstructures) {

    cat("\n ..",toupper(each_hem),"CORTEX ANALYSIS \n")
    verts_hem <- surf_list[[each_hem]]$vertices #verts_hem_old <- cifti_data[[each_hem]][[1]]$surf$vertices
    faces_hem <- surf_list[[each_hem]]$faces #faces_hem_old <- cifti_data[[each_hem]][[1]]$surf$faces

    #set up session list
    session_data <- vector('list', n_sess)
    names(session_data) <- session_names
    for(ss in 1:n_sess){
      sess <- list(BOLD = t(BOLD_list[[each_hem]][[ss]]), design=design[[ss]])
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
    mask_hem <- make_mask(session_data)
    BayesGLM_out <- BayesGLM(data = session_data,
                             beta_names = beta_names,
                             mesh = NULL,
                             vertices = verts_hem,
                             faces = faces_hem,
                             mask = mask_hem,
                             scale_BOLD = scale_BOLD,
                             scale_design = FALSE, # done above
                             Bayes = Bayes,
                             ar_order = ar_order,
                             ar_smooth = ar_smooth,
                             num.threads = num.threads,
                             return_INLA_result = return_INLA_result,
                             outfile = outfile_name,
                             verbose = verbose,
                             avg_sessions = avg_sessions,
                             trim_INLA = trim_INLA)

    BayesGLM_results[[each_hem]] <- BayesGLM_out[-grep("classical", names(BayesGLM_out))]
    classicalGLM_results[[each_hem]] <- BayesGLM_out$result_classical
    class(BayesGLM_results[[each_hem]]) <- 'BayesGLM'
    class(classicalGLM_results[[each_hem]]) <- 'classicalGLM'

    rm(BayesGLM_out); gc()
  }

  # update session info if averaged over sessions
  if(avg_sessions==TRUE){
    session_names <- 'session_avg'
    n_sess_orig <- n_sess
    n_sess <- 1
  } else {
    n_sess_orig <- NULL
  }

  ### CONSTRUCT BETA ESTIMATES AS CIFTI OBJECTS

  cat(' PUTTING RESULTS IN CIFTI FORMAT \n')

  classicalGLM_cifti <- BayesGLM_cifti <- vector('list', n_sess)
  names(classicalGLM_cifti) <- names(BayesGLM_cifti) <- session_names
  datL <- datR <- NULL
  for(ss in 1:n_sess){

    # CLASSICAL GLM
    if(do_left) datL <- classicalGLM_results$left[[ss]]$estimates[cortexL_mwall==1,]
    if(do_right) datR <- classicalGLM_results$right[[ss]]$estimates[cortexR_mwall==1,]
    classicalGLM_cifti[[ss]] <- as.xifti(
      cortexL = datL,
      cortexL_mwall = cortexL_mwall,
      cortexR = datR,
      cortexR_mwall = cortexR_mwall
    )
    classicalGLM_cifti[[ss]]$meta$cifti$names <- beta_names

    # BAYESIAN GLM
    if(do_Bayesian){
      if(do_left) datL <- BayesGLM_results$left$beta_estimates[[ss]][cortexL_mwall==1,]
      if(do_right) datR <- BayesGLM_results$right$beta_estimates[[ss]][cortexR_mwall==1,]
      BayesGLM_cifti[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = cortexL_mwall,
        cortexR = datR,
        cortexR_mwall = cortexR_mwall
      )
      BayesGLM_cifti[[ss]]$meta$cifti$names <- beta_names
    }
  }

  #TO DO: Combine hyperparameters across hemispheres, rename "beta" to "task" (and in BayesGLM)
  # Rename theta_posteriors to hyperpar_posteriors

  result <- list(betas_Bayesian = BayesGLM_cifti,
                 betas_classical = classicalGLM_cifti,
                 GLMs_Bayesian = list(cortexL = BayesGLM_results$left,
                                      cortexR = BayesGLM_results$right),
                 GLMs_classical = list(cortexL = classicalGLM_results$left,
                                       cortexR = classicalGLM_results$right),
                 session_names = session_names,
                 n_sess_orig = n_sess_orig,
                 beta_names = beta_names,
                 design = design) #task part of design matrix after centering/scaling but before nuisance regression and prewhitening

  cat('\n DONE! \n')

  class(result) <- "BayesGLM_cifti"
  return(result)
}
