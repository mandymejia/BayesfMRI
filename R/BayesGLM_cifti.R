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
#' @param nuisance (Optional) A TxJ matrix of nuisance signals
#'  (or list of such matrices, for multiple-session modeling).
#' @param dHRF Logical indicating whether the temporal derivative of each column
#'  in the design matrix should be added to \code{nuisance}. Default: \code{TRUE}.
#' @param hpf,DCT Add DCT bases to \code{nuisance} to apply a temporal
#'  high-pass filter to the data? Only one of these arguments should be provided.
#'  \code{hpf} should be the filter frequency; if it is provided, \code{TR}
#'  must be provided too. The number of DCT bases to include will be computed
#'  to yield a filter with as close a frequency to \code{hpf} as possible.
#'  Alternatively, \code{DCT} can be provided to directly specify the number
#'  of DCT bases to include.
#'
#'  Default: \code{DCT=4} (use four DCT bases for high-pass filtering; for
#'  typical \code{TR} this amounts to lower filter frequency than the
#'  approximately .01 Hz used in most studies.)
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @inheritParams Bayes_Param
#' @param ar_order (numeric) Controls prewhitening. If greater than zero, this
#'  should be a number indicating the order of the autoregressive model to use
#'  for prewhitening. If zero, do not prewhiten. Default: \code{6}.
#' @param ar_smooth FWHM parameter for smoothing. Remember that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} or \code{NULL}
#'  to not do any smoothing. Default: \code{5}.
#' @param aic Use the AIC to select AR model order between \code{0} and \code{ar_order}? Default: \code{FALSE}.
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
#' @param meanTol,varTol Tolerance for mean and variance of each data location. Locations which
#'  do not meet these thresholds are masked out of the analysis. Defaults: \code{1e-6}.
#' @param trim_INLA (logical) should the \code{INLA_result} objects within the
#'   result be trimmed to only what is necessary to use `id_activations()`? Default: `TRUE`.
#'
#' @return An object of class \code{"BayesGLM"}, a list containing...
#'
# @importFrom ciftiTools read_cifti resample_gifti as.xifti remove_xifti
#' @import ciftiTools
#' @importFrom matrixStats rowVars rowSums2 colVars
#' @importFrom parallel detectCores
#'
#' @export
BayesGLM_cifti <- function(
  cifti_fname,
  surfL_fname=NULL,
  surfR_fname=NULL,
  brainstructures=c('left','right'),
  design=NULL,
  onsets=NULL,
  TR=NULL,
  nuisance=NULL,
  dHRF=TRUE,
  hpf=NULL,
  DCT=if(is.null(hpf)) {4} else {NULL},
  scale_BOLD=c("auto", "mean", "sd", "none"),
  scale_design=TRUE,
  Bayes=TRUE,
  ar_order = 6,
  ar_smooth = 5,
  aic = FALSE,
  resamp_res=10000,
  num.threads=4,
  verbose=FALSE,
  outfile=NULL,
  return_INLA_result=FALSE,
  avg_sessions = TRUE,
  session_names=NULL,
  meanTol=1e-6, varTol=1e-6,
  trim_INLA = TRUE){

  # Check args ------------------------------------------
  is_xifti <- FALSE
  if (is.character(cifti_fname)) {
    NULL
  } else if (is.xifti(cifti_fname, messages=FALSE)) {
    is_xifti <- TRUE
    cifti_fname <- list(cifti_fname)
  } else if (is.list(cifti_fname)) {
    if (all(vapply(cifti_fname, is.character, FALSE)) && all(vapply(cifti_fname, length, 0)==1)) {
      cifti_fname <- as.character(cifti_fname)
    } else if (all(vapply(cifti_fname, is.xifti, messages=FALSE, FALSE))) {
      is_xifti <- TRUE
    }
  } else {
    stop('`cifti_fname` should be a character vector or list of `"xifti"` objects.')
  }

  # [TO DO]: If input is a `"xifti"`, infer `resamp_res`
  # or maybe just add surfaces to the `"xifti"` using `add_surf` and that will handle the
  # difference in resolution.

  # Rename and coerce to logical
  do_Bayesian <- as.logical(Bayes)
  if (do_Bayesian) { check_INLA(require_PARDISO=TRUE) }

  # Check nuisance arguments.
  stopifnot(is.logical(dHRF) && length(dHRF)==1)
  if (!is.null(DCT)) {
    stopifnot(is.numeric(DCT) && length(DCT)==1 && DCT>=0 && DCT==round(DCT))
    if (DCT==0) { DCT <- NULL }
  }
  if (!is.null(hpf)) {
    stopifnot(is.numeric(hpf) && length(hpf)==1 && hpf>=0)
    if (hpf==0) { hpf <- NULL }
  }

  # Check prewhitening arguments.
  if (is.null(ar_order)) ar_order <- 0
  ar_order <- as.numeric(ar_order)
  do_pw <- (ar_order > 0)
  if (is.null(ar_smooth)) ar_smooth <- 0
  ar_smooth <- as.numeric(ar_smooth)
  stopifnot(is.logical(aic) && length(aic)==1)

  avail_cores <- parallel::detectCores()
  num.threads <- min(num.threads, avail_cores)
  # if(avail_cores < 2) num.threads <- 1

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

  need_mesh <- (do_Bayesian || (do_pw && ar_smooth > 0))

  # Surfaces ---------------------------------------
  # if(do_left & is.null(surfL_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  # if(do_right & is.null(surfR_fname)) stop('surfL_fname must be provided if brainstructures includes "left"')
  surf_list <- list(left=NULL, right=NULL)
  if (need_mesh) {
    if (do_left) {
      if (is.null(surfL_fname)) {
        cat("Using `ciftiTools` default inflated surface for the left cortex.\n")
        surfL_fname <- ciftiTools.files()$surf["left"]
      }
      surf_list$left <- read_surf(surfL_fname, resamp_res=resamp_res)
    }
    if (do_right) {
      if (is.null(surfR_fname)) {
        cat("Using `ciftiTools` default inflated surface for the right cortex.\n")
        surfR_fname <- ciftiTools.files()$surf["right"]
      }
      surf_list$right <- read_surf(surfR_fname, resamp_res=resamp_res)
    }
  } else {
    surf_list <- list(
      left = list(vertices=NULL, faces=NULL),
      right = list(vertices=NULL, faces=NULL)
    )
  }

  # Sessions ---------------------------------------
  # Name sessions and check compatibility of multi-session arguments
  n_sess <- length(cifti_fname)
  if (n_sess == 1 & avg_sessions) avg_sessions <- FALSE
  if (n_sess==1) {
    if (is.null(session_names)) session_names <- 'single_session'
    if (!is.null(design)) design <- list(single_session = design)
    if (!is.null(onsets)) onsets <- list(single_session = onsets)
    if (!is.null(nuisance)) nuisance <- list(single_session = nuisance)
  } else {
    if (is.null(session_names)) session_names <- paste0('session', 1:n_sess)
    # if (length(session_names) == 1) { session_names <- paste0(session_names, 1:nsess) } # allow prefix?
    if (!is.null(design) && (length(design) != n_sess)) {
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`design` must be a list of length equal to the number of sessions ",
        " (or `NULL`, if onsets provided)."
      ))
    }
    if (!is.null(onsets) && (length(onsets) != n_sess)) {
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`onsets` must be a list of length equal to the number of sessions ",
        " (or `NULL`, if `design` provided)."
      ))
    }
    if (!is.null(nuisance) && (length(nuisance) != n_sess)) {
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`onsets` must be a list of length equal to the number of sessions ",
        " (or `NULL`)."
      ))
    }
  }
  if (length(session_names) != n_sess) {
    stop('If `session_names` is provided, it must be of the same length as `cifti_fname`')
  }

  # `design`, `onsets`, `beta_names` ---------------------------------------
  # Check `design` and `onsets`. Get `beta_names`.
  if (!xor(is.null(design), is.null(onsets))) { stop('`design` or `onsets` must be provided, but not both.') }
  if (!is.null(onsets) && is.null(TR)) { stop('Please provide `TR` if onsets provided') }
  if (!is.null(onsets)) {
    do_multisesh <- inherits(onsets[[1]], "list")
    if (!do_multisesh) { stopifnot(inherits(onsets, "matrix") || inherits(onsets, "data.frame")) }
    o1 <- if (do_multisesh) { onsets[[1]] } else { onsets }
    beta_names <- if (is.null(names(o1))) { paste0("beta", seq(length(o1))) } else { names(o1) }
    rm(o1)
  }
  if (!is.null(design)) {
    do_multisesh <- inherits(design, "list")
    if (!do_multisesh) { stopifnot(inherits(design, "matrix") || inherits(design, "data.frame")) }
    d1 <- if (do_multisesh) { design[[1]] } else { design }
    beta_names <- if (is.null(colnames(d1))) { paste0("beta", seq(ncol(d1))) } else { names(d1) }
    rm(d1)
  }

  # Scale
  if (isTRUE(scale_BOLD)) { cat("Setting `scale_BOLD <- 'auto'`"); scale_BOLD <- "auto" }
  if (isFALSE(scale_BOLD)) { cat("Setting `scale_BOLD <- 'none'`"); scale_BOLD <- "none" }
  scale_BOLD <- match.arg(scale_BOLD, c("auto", "mean", "sd", "none"))

  # Data setup ----------------------------------------------------------
  cat('\n SETTING UP DATA \n')
  ### For each session, separate the CIFTI data into left/right/sub and read in files
  BOLD_list <- list(left=NULL, right=NULL)
  mwallL <- mwallR <- NULL
  ntime <- vector("numeric", n_sess)

  for (ss in 1:n_sess) {
    if(n_sess > 1) cat(paste0(' .. reading in data for session ', ss,'\n'))

    if (is_xifti) {
      xii_ss <- cifti_fname[[ss]]
    } else {
      xii_ss <- read_cifti(
        cifti_fname[ss], brainstructures=brainstructures,
        resamp_res=resamp_res
      )
    }

    mwallL_ss <- xii_ss$meta$cortex$medial_wall_mask$left
    mwallR_ss <- xii_ss$meta$cortex$medial_wall_mask$right
    ntime[ss] <- ncol(xii_ss)

    # Get medial wall mask, or check that it matches.
    if (ss == 1) {
      if (do_left) { mwallL <- mwallL_ss }
      if (do_right) { mwallR <- mwallR_ss }
      # [TO DO] Check compatibility with `surf_list`
    } else {
      if (do_left) {
        stopifnot(length(mwallL) == length(mwallL_ss))
        stopifnot(all(mwallL == mwallL_ss))
      }
      if (do_right) {
        stopifnot(length(mwallR) == length(mwallR_ss))
        stopifnot(all(mwallR == mwallR_ss))
      }
    }

    # Grab BOLD data (input NAs in medial wall locations)
    if (do_left) {
      BOLD_list[["left"]][[ss]] <- ciftiTools::unmask_cortex(xii_ss$data$cortex_left, mwallL)
    }
    if (do_right) {
      BOLD_list[["right"]][[ss]] <- ciftiTools::unmask_cortex(xii_ss$data$cortex_right, mwallR)
    }
  }

  BOLD_list <- BOLD_list[!vapply(BOLD_list, is.null, FALSE)]
  if (need_mesh) { surf_list <- surf_list[!vapply(surf_list, is.null, FALSE)] }

  # Design and nuisance matrices -----------------------------------------------
  if (is.null(design)) {
    cat(" MAKING DESIGN MATRICES \n")
    design <- vector("list", n_sess)
    for (ss in seq(n_sess)) {
      design[[ss]] <- make_HRFs(onsets[[ss]], TR=TR, duration=ntime[ss])
    }
  }

  ### Check that design matrix names are consistent across sessions
  if (n_sess > 1) {
    tmp <- sapply(design, colnames)
    if(length(beta_names) == 1) {
      num_names <- length(unique(tmp))
      if (num_names > 1) stop('task names must match across sessions for multi-session modeling')
    } else {
      num_names <- apply(tmp, 1, function(x) length(unique(x))) #number of unique names per row
      if (max(num_names) > 1) stop('task names must match across sessions for multi-session modeling')
    }
  }

  if (scale_design) design <- sapply(design, scale_design_mat, simplify = F)
  if (!scale_design) design <- sapply(design, scale, scale = F, simplify = F) #center but do not scale

  ### ADD ADDITIONAL NUISANCE REGRESSORS
  DCTs <- vector("numeric", n_sess)
  for (ss in 1:n_sess) {
    # DCT highpass filter
    if (!is.null(hpf) || !is.null(DCT)) {
      # Get the num. of bases for this session.
      if (!is.null(hpf)) {
        DCTs[ss] <- round(dct_convert(ntime[ss], TR, f=hpf))
      } else {
        DCTs[ss] <- DCT
      }
      # Generate the bases and add them.
      DCTs[ss] <- dct_bases(ntime[ss], DCTs[ss])
      if (DCTs[ss] > 0) {
        if (!is.null(nuisance)) {
          nuisance[[ss]] <- cbind(nuisance[[ss]], DCTs[ss])
        } else {
          nuisance[[ss]] <- cbind(DCTs[ss])
        }
      }
    }
    # dHRF
    if (dHRF) {
      dHRF <- gradient(design[[ss]])
      if (!is.null(nuisance)) {
        nuisance[[ss]] <- cbind(nuisance[[ss]], dHRF)
      } else {
        nuisance[[ss]] <- dHRF
      }
    }
  }

  # Do GLM --------------------------------------------------------
  cat(' RUNNING MODELS \n')
  classicalGLM_results <- list(left = NULL, right = NULL)
  BayesGLM_results <- list(left = NULL, right = NULL)

  # >> Loop through brainstructures to complete the analyses on the different hemispheres ----
  for (each_hem in brainstructures) {

    cat("\n ..",toupper(each_hem),"CORTEX ANALYSIS \n")

    # set up session list
    session_data <- vector('list', n_sess)
    names(session_data) <- session_names
    for (ss in 1:n_sess) {
      session_data[[ss]] <- list(
        BOLD = t(BOLD_list[[each_hem]][[ss]]),
        design=design[[ss]]
      )
      if (!is.null(nuisance)) {
        session_data[[ss]]$nuisance <- nuisance[[ss]]
      }
    }

    # `outfile`
    if (!is.null(outfile)) {
      if (endsWith(outfile, ".rds")) {
        outfile_name <- gsub(".rds$", paste0("_",each_hem,".rds"), outfile)
      } else {
        outfile_name <- paste0(outfile, "_",each_hem,".rds")
      }
    } else {
      outfile_name <- NULL
    }

    BayesGLM_out <- BayesGLM(
      data = session_data,
      beta_names = beta_names,
      mesh = NULL,
      vertices = surf_list[[each_hem]]$vertices,
      faces = surf_list[[each_hem]]$faces,
      scale_BOLD = scale_BOLD,
      scale_design = FALSE, # done above
      Bayes = do_Bayesian,
      ar_order = ar_order,
      ar_smooth = ar_smooth,
      aic = aic,
      num.threads = num.threads,
      return_INLA_result = return_INLA_result,
      outfile = outfile_name,
      verbose = verbose,
      avg_sessions = avg_sessions,
      meanTol=meanTol,
      varTol=varTol,
      trim_INLA = trim_INLA
    )

    BayesGLM_results[[each_hem]] <- BayesGLM_out[-grep("classical", names(BayesGLM_out))]
    classicalGLM_results[[each_hem]] <- BayesGLM_out$result_classical
    class(BayesGLM_results[[each_hem]]) <- 'BayesGLM'
    class(classicalGLM_results[[each_hem]]) <- 'classicalGLM'

    rm(BayesGLM_out); gc()
  }

  # update session info if averaged over sessions
  if (avg_sessions==TRUE) {
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
    if(do_left) datL <- classicalGLM_results$left[[ss]]$estimates[mwallL==1,]
    if(do_right) datR <- classicalGLM_results$right[[ss]]$estimates[mwallR==1,]
    classicalGLM_cifti[[ss]] <- as.xifti(
      cortexL = datL,
      cortexL_mwall = mwallL,
      cortexR = datR,
      cortexR_mwall = mwallR
    )
    classicalGLM_cifti[[ss]]$meta$cifti$names <- beta_names

    # BAYESIAN GLM
    if(do_Bayesian){
      if(do_left) datL <- BayesGLM_results$left$beta_estimates[[ss]][mwallL==1,]
      if(do_right) datR <- BayesGLM_results$right$beta_estimates[[ss]][mwallR==1,]
      BayesGLM_cifti[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = mwallL,
        cortexR = datR,
        cortexR_mwall = mwallR
      )
      BayesGLM_cifti[[ss]]$meta$cifti$names <- beta_names
    }
  }

  #TO DO: Combine hyperparameters across hemispheres, rename "beta" to "task" (and in BayesGLM)
  # Rename theta_posteriors to hyperpar_posteriors

  result <- list(
    betas_Bayesian = BayesGLM_cifti,
    betas_classical = classicalGLM_cifti,
    GLMs_Bayesian = list(
      cortexL = BayesGLM_results$left,
      cortexR = BayesGLM_results$right
    ),
    GLMs_classical = list(
      cortexL = classicalGLM_results$left,
      cortexR = classicalGLM_results$right
    ),
    session_names = session_names,
    n_sess_orig = n_sess_orig,
    beta_names = beta_names,
    design = design
  ) #task part of design matrix after centering/scaling but before nuisance regression and prewhitening

  cat('\n DONE! \n')

  class(result) <- "BayesGLM_cifti"
  return(result)
}
