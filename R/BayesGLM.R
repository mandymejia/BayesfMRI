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
#' @param nuisance_include (Optional) Additional nuisance covariates to include.  Default is 'drift' (linear and quadratic drift terms) and 'dHRF' (temporal derivative of each column of design matrix).
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @param GLM_method One of \code{Bayesian}, \code{classical}, \code{EM}, or \code{all}, depending on the specific type of analysis that should be performed. Default method is \code{all}.
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
#' @param groups_df Data frame indicating the name and model group of each
#'   region.  See Details of the \code{BayesGLM_vol3D} function
#' @param tol Only used when \code{method} is set to 'EM'. The numeric tolerance
#'   used to determine convergence of the EM algorithm.
#'
#' @return An object of class \code{"BayesGLM"}, a list containing...
#'
# @importFrom ciftiTools read_cifti resample_gifti as.xifti remove_xifti
#' @import ciftiTools
#' @importFrom matrixStats rowVars rowSums2
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
                     GLM_method='all',
                     ar_order = 6,
                     ar_smooth = 5,
                     session_names=NULL,
                     resamp_res=10000,
                     num.threads=4,
                     verbose=FALSE,
                     outfile=NULL,
                     return_INLA_result=FALSE,
                     avg_sessions = TRUE,
                     trim_INLA = TRUE,
                     groups_df = NULL,
                     tol = NULL){

#  GLM_method = match.arg(GLM_method, c('both','Bayesian','classical'))
  GLM_method = match.arg(GLM_method, c('all','Bayesian','classical','EM'))

  do_Bayesian <- (GLM_method %in% c('all','Bayesian'))
  do_classical <- (GLM_method %in% c('all','classical'))
  do_EM <- (GLM_method %in% c('all','EM'))

  check_BayesGLM(require_PARDISO=(do_Bayesian|do_EM))

  avail_cores <- parallel::detectCores()
  if(avail_cores < 2) {
    num.threads <- 1
  }

  # Check that arguments are compatible
  brainstructures <- ciftiTools:::match_input(
    brainstructures, c("left","right","subcortical"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }
  # if(do_EM & "subcortical" %in% brainstructures)
    # stop("The EM algorithm does not yet work with the subcortical model.")
#   do_left <- ('left' %in% brainstructures)
#   do_right <- ('right' %in% brainstructures)
#   do_sub <- FALSE
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
      if(is.null(names(onsets)))
        beta_names <- paste0("beta",seq_len(length(onsets)))
      if(!is.null(names(onsets)))
        beta_names <- names(onsets)
    }
  }
  if(!is.null(design)) {
    if(any(class(design) == "list")) {
      if(is.null(colnames(design[[1]])))
        beta_names <- paste0("beta",seq_len(ncol(design[[1]])))
      if(!is.null(colnames(design[[1]])))
        beta_names <- colnames(design[[1]])
    }
    if("matrix" %in% class(design) | "data.frame" %in% class(design)) {
      if(is.null(colnames(design)))
        beta_names <- paste0("beta",seq_len(ncol(design)))
      if(!is.null(colnames(design)))
        beta_names <- colnames(design)
    }
  }

  # Check prewhitening arguments.
  if(is.null(ar_order)) ar_order <- 0
  ar_order <- as.numeric(ar_order)
  prewhiten <- (ar_order > 0)

  # if(do_left & is.null(surfL_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  # if(do_right & is.null(surfR_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  if('left' %in% brainstructures & is.null(surfL_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  if('right' %in% brainstructures & is.null(surfR_fname)) stop('surfR_fname must be provided if brainstructures includes "right"')

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

  cat('\nSETTING UP DATA\n')

  if(class(cifti_fname) == "character") {

    ### For each session, separate the CIFTI data into left/right/sub and read in files
    # if(do_left) cifti_left <- vector('list', n_sess)
    # if(do_right) cifti_right <- vector('list', n_sess)
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
    #     ntime <- ncol(cifti_right[[ss]])
    #   }
    #   #if(do_sub) { nifti_data[[ss]] <- cifti_ss$VOL; ntime <- ncol(cifti_ss$VOL) }
    #   #if(do_sub & ss==1) nifti_labels[[ss]] <- cifti_ss$LABELS

    #   if(make_design){
    #     cat(paste0('    Constructing design matrix for session ', ss, '\n'))
    #     design[[ss]] <- make_HRFs(onsets[[ss]], TR=TR, duration=ntime)
    #   }

    # }
    cat("READING IN CIFTI DATA \n")
    cifti_data <- sapply(brainstructures, function(br_str) {
      sapply(1:n_sess, function(ss){
        if(br_str %in% c("left","right")) {
          cifti_ss <- read_cifti(
            cifti_fname[ss],
            surfL_fname=surfL_fname, surfR_fname=surfR_fname,
            brainstructures=brainstructures,
            resamp_res=resamp_res
          )
          mwall_mask <- cifti_ss$meta$cortex$medial_wall_mask[[br_str]]
          cdata <- cifti_ss$data[[paste0("cortex_",br_str)]]
          cifti_out <- matrix(NA,
                              nrow = length(mwall_mask),
                              ncol=ncol(cdata))
          cifti_out[mwall_mask,] <- cdata
          return(list(
            data = cifti_out,
            surf = cifti_ss$surf[[paste0("cortex_", br_str)]],
            ntime = ncol(cdata),
            cifti = cifti_ss
          ))
        }
        if(br_str == "subcortical") {
          cifti_ss <- read_cifti(
            cifti_fname[ss],
            surfL_fname=NULL, surfR_fname=NULL,
            brainstructures=brainstructures,
            resamp_res=NULL
          )
          cdata <- cifti_ss$data$subcort
          return(list(
            data = cdata,
            surf = NULL,
            ntime = ncol(cdata),
            cifti = cifti_ss
          ))
        }
      }, simplify = F)
    }, simplify = F)
  }

  if(class(cifti_fname[[1]]) == "xifti") {
    cifti_data <- sapply(brainstructures, function(br_str) {
      if(br_str %in% c('left','right')) {
        hem_sess_cifti <- sapply(1:n_sess, function(ss){
          other_hem <- grep(br_str, c('left','right'), value = T, invert  = T)
          cifti_ss <- cifti_fname[[ss]]
          cifti_ss <- remove_xifti(cifti_ss, paste0("cortex_",other_hem))
          mwall_mask <- cifti_ss$meta$cortex$medial_wall_mask[[br_str]]
          cdata <- cifti_ss$data[[paste0("cortex_",br_str)]]
          cifti_out <- matrix(NA,
                              nrow = length(mwall_mask),
                              ncol=ncol(cdata))
          cifti_out[mwall_mask,] <- cdata
          return(list(
            data = cifti_out,
            surf = cifti_ss$surf[[paste0("cortex_", br_str)]],
            ntime = ncol(cdata),
            cifti = cifti_ss
          ))
        }, simplify = F)
        return(hem_sess_cifti)
      }
      if(br_str == "subcortical") {
        sub_sess_cifti <- sapply(1:n_sess, function(ss){
          cifti_ss <- cifti_fname[[ss]]
          cifti_ss <- remove_xifti(cifti_ss, paste0("cortex_",c('left','right')))
          cdata <- cifti_ss$data$subcort
          return(list(
            data = cdata,
            surf = NULL,
            ntime = ncol(cdata),
            cifti = cifti_ss
          ))
        }, simplify = F)
        return(sub_sess_cifti)
      }
    }, simplify = F)
  }

  ntimes <- sapply(cifti_data, function(h) sapply(h, `[[`, i = 3), simplify = F)

  if(is.null(design)) {
    cat("MAKING DESIGN MATRICES \n")
    design <- mapply(make_HRFs, onsets, TR = TR, duration = ntimes[[1]], SIMPLIFY = F)
  }

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

  # classicalGLM_left <- classicalGLM_right <- classicalGLM_vol <- NULL
  # BayesGLM_left <- BayesGLM_right <- BayesGLM_vol <- NULL
  classicalGLM_results <- list(left=NULL,right=NULL,vol=NULL)
  BayesGLM_results <- list(left=NULL,right=NULL,vol=NULL)
  GLMEM_results <- list(left=NULL, right=NULL,vol=NULL)
  pw_data <- list(left=NULL,right=NULL,vol=NULL)

  ### FORMAT DESIGN MATRIX
  # for(ss in 1:n_sess){
  #   if(scale_design){
  #     design[[ss]] <- scale_design_mat(design[[ss]])
  #   } else {
  #     design[[ss]] <- scale(design[[ss]], scale=FALSE) #center design matrix to eliminate baseline
  #   }
  # }
  if(scale_design) design <- sapply(design, scale_design_mat, simplify = F)
  if(!scale_design) design <- sapply(design, scale, scale=F, simplify = F)

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
  for(br_str in brainstructures) {
    # >>> Cortex models ----
    if(br_str %in% c("left","right")) {
      cat("\n ...",toupper(br_str),"CORTEX \n")
      verts <- cifti_data[[br_str]][[1]]$surf$vertices
      faces <- cifti_data[[br_str]][[1]]$surf$faces

      #set up session list
      session_data <- vector('list', n_sess)
      names(session_data) <- session_names
      for(ss in 1:n_sess){
        sess <- list(BOLD = t(cifti_data[[br_str]][[ss]]$data), design=design[[ss]])
        if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
        session_data[[ss]] <- sess
      }

      scale_BOLD_hem <- scale_BOLD

      if(prewhiten) {
        pw_data[[br_str]] <-
          prewhiten_cifti(
            data = session_data,
            mask = NULL,
            scale_BOLD = TRUE,
            scale_design = TRUE,
            ar_order = ar_order,
            ar_smooth = ar_smooth,
            cifti_data = cifti_data[[br_str]][[1]]$cifti,
            brainstructure = br_str,
            num.threads = num.threads
          )
        session_data <- pw_data[[br_str]]$data
        scale_BOLD_hem <- FALSE # Done above
      }

      if(!is.null(outfile)) {
        if (endsWith(outfile, ".rds")) {
          outfile_name <- gsub(".rds$", paste0("_",br_str,".rds"), outfile)
        } else {
          outfile_name <- paste0(outfile, "_",br_str,".rds")
        }
      } else {
        outfile_name <- NULL
      }

      other_hem <- grep(br_str, c("left","right"), value = T, invert = T)
      if(do_classical) {
        # Adding the try() function here so that if one model fails, the others
        # can still run
        start_time <- proc.time()[3]
        classicalGLM_results[[br_str]] <- try(classicalGLM(data=session_data,
                                                       scale_BOLD=scale_BOLD_hem,
                                                       scale_design = FALSE)) # done above
        classicalGLM_results[[br_str]]$total_time <- proc.time()[3] - start_time
      }
      if(do_Bayesian) {
        start_time <- proc.time()[3]
        BayesGLM_results[[br_str]] <- try(BayesGLM(data = session_data,
                                               beta_names = beta_names,
                                               vertices = verts,
                                               faces = faces,
                                               scale_BOLD = scale_BOLD_hem,
                                               scale_design = FALSE, # done above
                                               num.threads = num.threads,
                                               return_INLA_result = return_INLA_result,
                                               outfile = outfile_name,
                                               verbose = verbose,
                                               avg_sessions = avg_sessions,
                                               trim_INLA = trim_INLA))
        BayesGLM_results[[br_str]]$total_time <- proc.time()[3] - start_time
      }
      if(do_EM) {
        start_time <- proc.time()[3]
        GLMEM_results[[br_str]] <- try(BayesGLMEM(data = session_data,
                                              beta_names = beta_names,
                                              vertices = verts,
                                              faces = faces,
                                              scale_BOLD = scale_BOLD_hem,
                                              scale_design = FALSE,
                                              EM_method = "separate",
                                              use_SQUAREM = TRUE,
                                              tol = tol,
                                              num.threads = num.threads,
                                              outfile = outfile_name,
                                              verbose = verbose))
        GLMEM_results[[br_str]]$total_time <- proc.time()[3] - start_time
      }
    }
    # >>> Subcortical model ----
    if(br_str == "subcortical") {
      cat("\n ... SUBCORTEX \n")
      locs <- which(cifti_data[[br_str]][[1]]$cifti$meta$subcort$mask, arr.ind = T)
      labs <- cifti_data[[br_str]][[1]]$cifti$meta$subcort$labels

      #set up session list
      session_data <- vector('list', n_sess)
      names(session_data) <- session_names
      for(ss in 1:n_sess){
        sess <- list(BOLD = t(cifti_data[[br_str]][[ss]]$data), design=design[[ss]])
        if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
        session_data[[ss]] <- sess
      }

      scale_BOLD_sub <- scale_BOLD

      # THIS STILL NEEDS TO BE TESTED
      if(prewhiten) {
        pw_data[[br_str]] <-
          prewhiten_cifti(
            data = session_data,
            mask = NULL,
            scale_BOLD = TRUE,
            scale_design = TRUE,
            ar_order = ar_order,
            ar_smooth = ar_smooth,
            cifti_data = cifti_data[[br_str]][[1]]$cifti,
            brainstructure = br_str,
            num.threads = num.threads
          )
        session_data <- pw_data[[br_str]]$data
        scale_BOLD_sub <- FALSE # Done above
      }

      if(!is.null(outfile)) {
        if (endsWith(outfile, ".rds")) {
          outfile_name <- gsub(".rds$", paste0("_",br_str,".rds"), outfile)
        } else {
          outfile_name <- paste0(outfile, "_",br_str,".rds")
        }
      } else {
        outfile_name <- NULL
      }

      if(do_classical) {
        start_time <- proc.time()[3]
        classicalGLM_results$vol <- try(classicalGLM(data=session_data,
                                                       scale_BOLD=scale_BOLD_sub,
                                                       scale_design = FALSE)) # done above
        classicalGLM_results$vol$total_time <- proc.time()[3] - start_time
      }
      if(do_Bayesian) {
        start_time <- proc.time()[3]
        BayesGLM_results$vol <- try(BayesGLM_vol3D(
          data = session_data,
          locations = locs,
          labels = labs,
          groups_df = groups_df,
          scale_BOLD=TRUE,
          scale_design = FALSE, # Done above
          return_INLA_result=return_INLA_result,
          outfile = outfile_name,
          num.threads = num.threads,
          verbose=verbose
        ))
        BayesGLM_results$vol$total_time <- proc.time()[3] - start_time
      }
      if(do_EM) {
        start_time <- proc.time()[3]
        GLMEM_results$vol <- try(BayesGLMEM_vol3D(
          data = session_data,
          beta_names = NULL,
          locations = locs,
          labels = labs,
          EM_method = "separate",
          use_SQUAREM = TRUE,
          groups_df = groups_df,
          scale_BOLD=TRUE,
          scale_design = FALSE, # Done above
          tol = tol,
          outfile = outfile_name,
          num.threads = num.threads,
          verbose=verbose,
          avg_sessions = TRUE
        ))
        GLMEM_results$vol$total_time <- proc.time()[3] - start_time
      }
    }
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

  classicalGLM_cifti <- BayesGLM_cifti <- GLMEM_cifti <- vector('list', n_sess)
  names(classicalGLM_cifti) <- names(BayesGLM_cifti) <- names(GLMEM_cifti) <- session_names
  datL <- datR <- dat_sub <- cortexL_mwall <- cortexR_mwall <- sub_mask <- out_labs <- NULL
  do_left <- 'left' %in% brainstructures
  do_right <- 'right' %in% brainstructures
  do_sub <- "subcortical" %in% brainstructures
  if(do_left) cortexL_mwall <- as.numeric(cifti_data$left[[1]]$cifti$meta$cortex$medial_wall_mask$left)
  if(do_right) cortexR_mwall <- as.numeric(cifti_data$right[[1]]$cifti$meta$cortex$medial_wall_mask$right)
  if(do_sub) {
    sub_mask <- cifti_data$subcortical[[1]]$cifti$meta$subcort$mask
    out_labs <- labs
  }
  for(ss in 1:n_sess){
    if(do_classical){
      if(do_left & class(classicalGLM_results$left) != "try-error") {
        datL <- classicalGLM_results$left[[ss]]$estimates[cortexL_mwall==1,]
      }
      if(do_right & class(classicalGLM_results$right) != "try-error") {
        datR <- classicalGLM_results$right[[ss]]$estimates[cortexR_mwall==1,]
      }
      if(do_sub & class(classicalGLM_results$vol) != "try-error") {
        dat_sub <- classicalGLM_results$vol[[ss]]$estimates
      }
      classicalGLM_cifti[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = cortexL_mwall,
        cortexR = datR,
        cortexR_mwall = cortexR_mwall,
        subcortVol = dat_sub,
        subcortLabs = out_labs,
        subcortMask = sub_mask
        #subcortVol = classicalGLM_vol$single_session,
        #mask = mask,
        #subcortLab = nifti_labels
      )

      classicalGLM_cifti[[ss]]$meta$cifti$names <- beta_names
    }
    if(do_Bayesian){
      if(do_left & class(BayesGLM_results$left) != "try-error")
        datL <- BayesGLM_results$left$beta_estimates[[ss]][cortexL_mwall==1,]
      if(do_right & class(BayesGLM_results$right) != "try-error")
        datR <- BayesGLM_results$right$beta_estimates[[ss]][cortexR_mwall==1,]
      if(do_sub & class(BayesGLM_results$vol) != "try-error")
        dat_sub <- BayesGLM_results$vol$beta_estimates[[ss]]
      BayesGLM_cifti[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = cortexL_mwall,
        cortexR = datR,
        cortexR_mwall = cortexR_mwall,
        subcortVol = dat_sub,
        subcortLabs = out_labs,
        subcortMask = sub_mask
        #subcortVol = BayesGLM_vol$single_session,
        #mask = mask,
        #subcortLab = nifti_labels
      )
      BayesGLM_cifti[[ss]]$meta$cifti$names <- beta_names
    }
    if(do_EM) {
      if(do_left & class(GLMEM_results$left) != "try-error")
        datL <- GLMEM_results$left$beta_estimates[[ss]]
      if(do_right & class(GLMEM_results$right) != "try-error")
        datR <- GLMEM_results$right$beta_estimates[[ss]]
      if(do_sub & class(GLMEM_results$vol) != "try-error")
        dat_sub <- GLMEM_results$vol$beta_estimates[[ss]]
      GLMEM_cifti[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = cortexL_mwall,
        cortexR = datR,
        cortexR_mwall = cortexR_mwall,
        subcortVol = dat_sub,
        subcortLabs = out_labs,
        subcortMask = sub_mask
        #subcortVol = BayesGLM_vol$single_session,
        #mask = mask,
        #subcortLab = nifti_labels
      )
      GLMEM_cifti[[ss]]$meta$cifti$names <- beta_names
    }
  }

  if(avg_sessions) {
    if(do_classical) {
      # Need to include the missing locations for medial walls within the
      # linear combinations or the xifti object will be mis-mapped.
      if(do_left & class(classicalGLM_results$left) != "try-error") {
        datL <- matrix(NA, sum(cortexL_mwall),length(beta_names))
        maskL <- classicalGLM_results$left$avg$mask[cortexL_mwall == 1]
        datL[maskL,] <- classicalGLM_results$left$avg$estimates[cortexL_mwall == 1,]
      }
      if(do_right & class(classicalGLM_results$right) != "try-error") {
        datR <- matrix(NA, sum(cortexR_mwall),length(beta_names))
        maskR <- classicalGLM_results$right$avg$mask[cortexR_mwall == 1]
        datR[maskR,] <- classicalGLM_results$right$avg$estimates[cortexR_mwall == 1,]
      }
      if(do_sub & class(classicalGLM_results$vol) != "try-error") {
        dat_sub <- classicalGLM_results$vol$avg$estimates
      }
      # Adding the averages to the front of the BayesGLM_cifti object
      classicalGLM_cifti <- c(
        list(avg = as.xifti(
          cortexL = datL,
          cortexL_mwall = cortexL_mwall,
          cortexR = datR,
          cortexR_mwall = cortexR_mwall,
          subcortVol = dat_sub,
          subcortLabs = out_labs,
          subcortMask = sub_mask
        )),
        classicalGLM_cifti
      )
      classicalGLM_cifti[[1]]$meta$cifti$names <- beta_names
    }

    if (do_Bayesian) {
      # Need to include the missing locations for medial walls within the
      # linear combinations or the xifti object will be mis-mapped.
      if(do_left & class(BayesGLM_results$left) != "try-error") {
        datL <- matrix(NA, sum(cortexL_mwall),length(beta_names))
        maskL <- BayesGLM_results$left$mask[cortexL_mwall == 1]
        datL[maskL,] <- BayesGLM_results$left$avg_beta_estimates
      }
      if(do_right & class(BayesGLM_results$right) != "try-error") {
        datR <- matrix(NA, sum(cortexR_mwall),length(beta_names))
        maskR <- BayesGLM_results$right$mask[cortexR_mwall == 1]
        datR[maskR,] <- BayesGLM_results$right$avg_beta_estimates
      }
      if(do_sub & class(BayesGLM_results$vol) != "try-error") {
        dat_sub <- BayesGLM_results$vol$avg_beta_estimates
      }
      # Adding the averages to the front of the BayesGLM_cifti object
      BayesGLM_cifti <- c(
        list(avg = as.xifti(
          cortexL = datL,
          cortexL_mwall = cortexL_mwall,
          cortexR = datR,
          cortexR_mwall = cortexR_mwall,
          subcortVol = dat_sub,
          subcortLabs = out_labs,
          subcortMask = sub_mask
        )),
        BayesGLM_cifti
      )
      BayesGLM_cifti[[1]]$meta$cifti$names <- beta_names
    }

    if (do_EM) {
      # Need to include the missing locations for medial walls within the
      # linear combinations or the xifti object will be mis-mapped.
      if(do_left & class(GLMEM_results$left) != "try-error") {
        datL <- matrix(NA, sum(cortexL_mwall),length(beta_names))
        maskL <- GLMEM_results$left$mask[cortexL_mwall == 1]
        datL[maskL,] <- GLMEM_results$left$avg_beta_estimates
      }
      if(do_right & class(GLMEM_results$right) != "try-error") {
        datR <- matrix(NA, sum(cortexR_mwall),length(beta_names))
        maskR <- GLMEM_results$right$mask[cortexR_mwall == 1]
        datR[maskR,] <- GLMEM_results$right$avg_beta_estimates
      }
      if(do_sub & class(GLMEM_results$vol) != "try-error") {
        dat_sub <- GLMEM_results$vol$avg_beta_estimates
      }
      # Adding the averages to the front of the BayesGLM_cifti object
      GLMEM_cifti <- c(
        list(avg = as.xifti(
          cortexL = datL,
          cortexL_mwall = cortexL_mwall,
          cortexR = datR,
          cortexR_mwall = cortexR_mwall,
          subcortVol = dat_sub,
          subcortLabs = out_labs,
          subcortMask = sub_mask
        )),
        GLMEM_cifti
      )
      GLMEM_cifti[[1]]$meta$cifti$names <- beta_names
    }
  }


  prewhitening_info <- list()
  if(prewhiten) {
    # I take out the first element here because the first argument
    # is a copy of the data, which are already provided.
    if(do_left) prewhitening_info$left <- pw_data$left[-1]
    if(do_right) prewhitening_info$right <- pw_data$right[-1]
  }

  result <- list(session_names = session_names,
                 beta_names = beta_names,
                 betas_Bayesian = BayesGLM_cifti,
                 betas_classical = classicalGLM_cifti,
                 betas_EM = GLMEM_cifti,
                 GLMs_Bayesian = list(cortexL = BayesGLM_results$left,
                                      cortexR = BayesGLM_results$right,
                                      subcortical = BayesGLM_results$vol),
                 GLMs_classical = list(cortexL = classicalGLM_results$left,
                                       cortexR = classicalGLM_results$right,
                                       subcortical = classicalGLM_results$vol),
                 GLMs_EM = list(cortexL = GLMEM_results$left,
                                cortexR = GLMEM_results$right,
                                subcortical = GLMEM_results$vol),
                 prewhitening_info = prewhitening_info,
                 design = design)

  cat('\n DONE! \n')

  class(result) <- "BayesGLM_cifti"
  return(result)
}



#' BayesGLM
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
#' @inheritParams num.threads_Param
#' @inheritParams return_INLA_result_Param_TRUE
#' @param outfile File name where results will be written (for use by
#'  \code{BayesGLM2}).
#' @inheritParams verbose_Param_inla
#' @inheritParams contrasts_Param_inla
#' @inheritParams avg_sessions_Param
#' @param trim_INLA (logical) should the \code{INLA_result} objects within the
#'   result be trimmed to only what is necessary to use `id_activations()`? Default: `TRUE`.
#'
#' @return A list containing...
#'
#' @importFrom INLA inla.spde2.matern inla.pardiso.check inla.setOption inla.make.lincombs
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#'
#' @export
BayesGLM <- function(
  data,
  beta_names = NULL,
  vertices = NULL,
  faces = NULL,
  mesh = NULL,
  mask = NULL,
  scale_BOLD=TRUE,
  scale_design = TRUE,
  num.threads=4,
  return_INLA_result=TRUE,
  outfile = NULL,
  verbose=FALSE,
  contrasts = NULL,
  avg_sessions = TRUE,
  trim_INLA = TRUE){

  #check whether data is a list OR a session (for single-session analysis)
  #check whether each element of data is a session (use is.session)
  # V = number of data locations
  # T = length of time series for each session (vector)
  # K = number of unique tasks in all sessions

  check_BayesGLM(require_PARDISO=TRUE)

  #check that only mesh OR vertices+faces supplied
  has_mesh <- !is.null(mesh)
  has_verts_faces <- !is.null(vertices) & !is.null(faces)
  has_howmany <- has_mesh + has_verts_faces
  if(has_howmany != 1) stop('Must supply EITHER mesh OR vertices and faces.')

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data)
  n_sess <- length(session_names)
  if(n_sess == 1 & avg_sessions) avg_sessions <- FALSE

  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  # V <- ncol(data[[1]]$BOLD) #number of data locations
  # K <- ncol(data[[1]]$design) #number of tasks
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


  #Check that data locations same across sessions
  if(n_sess > 1){
    is_missing_all <- sapply(data, function(x) is.na(x$BOLD[1,]))
    tmp <- is_missing_all - is_missing
    if(max(abs(tmp))>0) stop('Missing (NA) data locations in BOLD must be consistent across sessions, but they are not.')
  }

  # for(s in 1:n_sess){
  #   if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
  #   if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
  #   if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  # }
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

  #remove zero var locations from mask
  if(sum(zero_var) > 0){
    if(is.null(mask)) num_flat <- sum(zero_var) else num_flat <- sum(zero_var[mask==1])
    if(num_flat > 1) warning(paste0('I detected ', num_flat, ' vertices that are flat (zero variance), NA or NaN in at least one session. Removing these from analysis. See mask returned with function output.'))
    if(num_flat == 1) warning(paste0('I detected 1 vertex that is flat (zero variance), NA or NaN in at least one session. Removing it from analysis. See mask returned with function output.'))
    mask_orig <- mask
    if(!is.null(mask)) mask[zero_var==TRUE] <- 0
    if(is.null(mask)) mask <- !zero_var
  } else {
    mask_orig <- NULL
  }

  #apply mask to mesh
  if(!is.null(mask)) {
    mask <- as.logical(mask)
    mesh_orig <- mesh #for later plotting
    mesh <- excursions::submesh.mesh(mask, mesh)
    mask <- !is.na(mesh$idx$loc) #update mask (sometimes vertices not excluded by mask will be excluded in mesh)
    mesh$idx$loc <- mesh$idx$loc[mask]
    for(s in 1:n_sess){
      data[[s]]$BOLD <- data[[s]]$BOLD[,mask]
    }
    V <- sum(mask)
    #zero_var <- zero_var[mask]
  } else {
    mesh_orig <- NULL
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
  replicates_list <- organize_replicates(n_sess=n_sess, beta_names=beta_names, mesh=mesh)
  betas <- replicates_list$betas
  repls <- replicates_list$repls

  #organize the formula and data objects
  #formula <- make_formula(beta_names = names(betas), repl_names = names(repls), hyper_initial = c(-2,2))
  #formula <- as.formula(formula)

  repl_names <- names(repls)
  hyper_initial <- c(-2,2)
  hyper_initial <- rep(list(hyper_initial), K)
  hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

  formula_vec <- paste0('f(',beta_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
  formula_vec <- c('y ~ -1', formula_vec)
  formula_str <- paste(formula_vec, collapse=' + ')
  formula <- as.formula(formula_str, env = globalenv())

  model_data <- make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)

  if(n_sess > 1 & avg_sessions) {
    cat("Set linear combinations for averages across sessions\n")
    diag_coefs <- Diagonal(n = V, x = 1/n_sess)
    session_coefs <- Matrix::bdiag(rep(list(diag_coefs),length(beta_names))) # Just finished this line
    full_coefs <- Reduce(cbind,rep(list(session_coefs),n_sess))
    design_pred <- rbind(model_data$X, full_coefs)
    response_pred <- rep(NA,dim(full_coefs)[1])
    model_data$y <- c(model_data$y,response_pred)
    #estimate model using INLA
    cat('\n ...... estimating model with INLA')
    system.time(
      INLA_result <- estimate_model(
        formula=formula, data=model_data, A=design_pred, spde, prec_initial=1,
        num.threads=num.threads, verbose=verbose, contrasts = contrasts
      )
    )
    cat('\n ...... model estimation completed')
  } else {
    #estimate model using INLA
    cat('\n ...... estimating model with INLA')
    system.time(
      INLA_result <- estimate_model(
        formula=formula, data=model_data, A=model_data$X, spde, prec_initial=1,
        num.threads=num.threads, verbose=verbose, contrasts = contrasts
      )
    )
    cat('\n ...... model estimation completed')
  }

  #extract useful stuff from INLA model result
  beta_estimates <- extract_estimates(object=INLA_result, session_names=session_names, mask=mask) #posterior means of latent task field
  theta_posteriors <- get_posterior_densities(object=INLA_result, spde, beta_names) #hyperparameter posterior densities
  # The mean of the mean beta estimates across sessions
  if(n_sess > 1 & avg_sessions) {
    pred_idx <- which(is.na(model_data$y))
    INLA_result$misc$configs$config[[1]]$pred_idx <- pred_idx
    avg_beta_means <- INLA_result$summary.linear.predictor$mean[pred_idx]
    avg_beta_estimates <- sapply(seq(K), function(k) {
      bbeta_out <- avg_beta_means[(seq(V) + (k-1)*V)]
      return(bbeta_out)
    }, simplify = T)
    names(avg_beta_estimates) <- beta_names
  } else {
    avg_beta_estimates <- NULL
  }

  #extract stuff needed for group analysis
  mu.theta <- INLA_result$misc$theta.mode
  Q.theta <- solve(INLA_result$misc$cov.intern)

  #construct object to be returned
  if(return_INLA_result){
    if(trim_INLA) INLA_result <- trim_INLA_result(INLA_result)
    result <- list(INLA_result = INLA_result,
                   mesh = mesh,
                   mesh_orig = mesh_orig,
                   mask = mask,
                   mask_orig = mask_orig,
                   design = design,
                   session_names = session_names,
                   beta_names = beta_names,
                   beta_estimates = beta_estimates,
                   avg_beta_estimates = avg_beta_estimates,
                   theta_posteriors = theta_posteriors,
                   mu.theta = mu.theta, #for joint group model
                   Q.theta = Q.theta, #for joint group model
                   y = y_all, #for joint group model
                   X = X_all_list, #for joint group model
                   #model_data, #temporary
                   #formula, #temporary
                   call = match.call())
  } else {
    result <- list(INLA_result = NULL,
                   mesh = mesh,
                   mesh_orig = mesh_orig,
                   mask = mask,
                   mask_orig = mask_orig,
                   design = design,
                   session_names = session_names,
                   beta_names = beta_names,
                   beta_estimates = beta_estimates,
                   avg_beta_estimates = avg_beta_estimates,
                   theta_posteriors = theta_posteriors,
                   mu.theta = mu.theta, #for joint group model
                   Q.theta = Q.theta, #for joint group model
                   y = y_all, #for joint group model
                   X = X_all_list, #for joint group model
                   #model_data, #temporary
                   #formula, #temporary
                   call = match.call())
  }


  class(result) <- "BayesGLM"

  if(!is.null(outfile)){
    saveRDS(result, file=outfile)
  }

  return(result)
}
