#' BayesGLM for CIFTI
#'
#' Performs spatial Bayesian GLM on the cortical surface for fMRI task activation
#'
#' @section INLA latent fields limit:
#'  INLA computation times increase greatly when the number of columns in the
#'  design matrix exceeds five. So if there are more than five tasks, or
#'  three or more tasks each with its temporal derivative being modeled as a
#'  task, \code{BayesGLM} will raise a warning. In cases like the latter, we
#'  recommend modeling the temporal derivatives as nuisance signals using the
#'  \code{nuisance} argument, rather than modeling them as tasks.
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @section Connectome Workbench Requirement:
#'  This function uses a system wrapper for the 'wb_command' executable. The
#'  user must first download and install the Connectome Workbench, available
#'  from https://www.humanconnectome.org/software/get-connectome-workbench .
#'
#' @param cifti_fname fMRI timeseries data in CIFTI format ("*.dtseries.nii").
#'  For single-session analysis this can be a file path to a CIFTI file or a
#'  \code{"xifti"} object from the \code{ciftiTools} package. For multi-session
#'  analysis this can be a vector of file paths or a list of \code{"xifti"}
#'  objects.
#' @param surfL_fname Left cortex surface geometry in GIFTI format
#'  ("*.surf.gii"). This can be a file path to a GIFTI file or a \code{"surf"}
#'  object from the \code{ciftiTools} package. This argument is only used if
#'  \code{brainstructures} includes \code{"left"} and \code{Bayes==TRUE}. If
#'  it's not provided, the HCP group-average inflated surface included in the
#'  \code{ciftiTools} package will be used.
#' @param surfR_fname Right cortex surface geometry in GIFTI format
#'  ("*.surf.gii"). This can be a file path to a GIFTI file or a \code{"surf"}
#'  object from the \code{ciftiTools} package. This argument is only used if
#'  \code{brainstructures} includes \code{"right"} and \code{Bayes==TRUE}. If
#'  it's not provided, the HCP group-average inflated surface included in the
#'  \code{ciftiTools} package will be used.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to analyze: \code{"left"} (left cortical surface) and/or \code{"right"}
#'  (right cortical surface). Default: \code{c("left","right")} (both
#'  hemispheres). Note that the subcortical models have not yet been implemented.
#' @param design,onsets,TR Either provide \code{design} directly, or provide
#'  both \code{onsets} and \code{TR} from which the design matrix or matrices
#'  will be constructed.
#'
#'   \code{design} is a \eqn{T \times K} task design matrix. Each column
#'   represents the expected BOLD response due to each task, a convolution of
#'   the hemodynamic response function (HRF) and the task stimulus. Note that
#'   the scale of the regressors will affect the scale and interpretation of the
#'   beta coefficients, so imposing a proper scale is recommended; see the
#'   \code{scale_design} argument, which by default is \code{TRUE}. Task names
#'   should be the column names, if not provided by the \code{task_names}
#'   argument. For multi-session modeling, this argument should be a list of
#'   such matrices. To model HRF derivatives, calculate the derivatives of the
#'   task columns beforehand (see the helper function \code{\link{cderiv}} which
#'   computes the discrete central derivative) and either add them to
#'   \code{design} to model them as tasks, or \code{nuisance} to model them as
#'   nuisance signals; it's recommended to then drop the first and last
#'   timepoints because the discrete central derivative doesn't exist at the
#'   time series boundaries. Do note that INLA computation times increase
#'   greatly if the design matrix has more than five columns, so it might be
#'   required to add these derivatives to \code{nuisance} rather than
#'   \code{design}.
#'
#'   \code{onsets} is an \eqn{L}-length list in which the name of each element is
#'   the name of the corresponding task, and the value of each element is a
#'   matrix of onsets (first column) and durations (second column) for each
#'   stimuli (each row) of the corresponding task. The units of both columns
#'   is seconds. For multi-session modeling, this argument should be a list of
#'   such lists. To model HRF derivatives, use the arguments \code{dHRF} and
#'   \code{dHRF_as}. If \code{dHRF==0} or \code{dHRF_as=="nuisance"}, the total
#'   number of columns in the design matrix, \eqn{K}, will equal \eqn{L}.
#'   If \code{dHRF_as=="task"}, \eqn{K} will equal \eqn{L} times \code{dHRF+1}.
#'
#'   \code{TR} is the temporal resolution of the data, in seconds.
#'
#' @param nuisance (Optional) A \eqn{T \times J} matrix of nuisance signals.
#'  These are regressed from the fMRI data and the design matrix prior to the
#'  GLM computation. For multi-session modeling, this argument should be a list
#'  of such matrices.
#' @param dHRF,dHRF_as Only applicable if \code{onsets} and \code{TR} are
#'  provided. These arguments enable the modeling of HRF derivatives.
#'
#'  Set \code{dHRF} to \code{1} to model the temporal derivatives of each task,
#'  \code{2} to add the second derivatives too, or \code{0} to not model the
#'  derivatives. Default: \code{1}.
#'
#'  If \code{dHRF > 0}, \code{dHRF_as} controls whether the derivatives are
#'  modeled as \code{"nuisance"} signals to regress out, \code{"tasks"}, or
#'  \code{"auto"} (default) to treat them as tasks unless the total number of
#'  columns in the design matrix would exceed five.
#'
#' @param hpf,DCT Add DCT bases to \code{nuisance} to apply a temporal
#'  high-pass filter to the data? Only one of these arguments should be provided.
#'  \code{hpf} should be the filter frequency; if it is provided, \code{TR}
#'  must be provided too. The number of DCT bases to include will be computed
#'  to yield a filter with as close a frequency to \code{hpf} as possible.
#'  Alternatively, \code{DCT} can be provided to directly specify the number
#'  of DCT bases to include.
#'
#'  Default: \code{DCT=4}. For typical \code{TR}, four DCT bases amounts to a
#'  lower frequency cutoff than the approximately .01 Hz used in most studies.
#'  We selected this default to err on the side of retaining more low-frequency
#'  information, but we recommend setting these arguments to values most
#'  appropriate for the data analysis at hand.
#'
#'  Using at least two DCT bases is as sufficient as using linear and quadratic
#'  drift terms in the design matrix. So if DCT detrending is being used, there
#'  is no need to add linear and quadratic drift terms to \code{nuisance}.
#' @param resamp_res The number of vertices to which each cortical surface
#'  should be resampled, or \code{NULL} to not resample. For computational
#'  feasibility, a value of \code{10000} or lower is recommended.
#' @inheritParams task_names_Param
#' @inheritParams session_names_Param
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @inheritParams Bayes_Param
# @inheritParams EM_Param
#' @inheritParams ar_order_Param
#' @inheritParams ar_smooth_Param
#' @inheritParams aic_Param
#' @inheritParams num.threads_Param
#' @inheritParams return_INLA_Param
#' @inheritParams verbose_Param
#' @inheritParams combine_sessions_Param
#' @param meanTol,varTol Tolerance for mean and variance of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Default: \code{1e-6} for both.
# @inheritParams emTol_Param
#' @inheritParams trim_INLA_Param
#'
#' @return An object of class \code{"BayesGLM_cifti"}: a list with elements
#'  \describe{
#'    \item{betas_Bayesian}{The task coefficients for the Bayesian model.}
#'    \item{betas_classical}{The task coefficients for the classical model.}
#'    \item{GLMs_Bayesian}{The entire list of GLM results, except for parameters estimated for the classical model.}
#'    \item{GLMs_classical}{Parameters estimated for the classical model from the GLM.}
#'    \item{session_names}{The names of the sessions.}
#'    \item{n_sess_orig}{The number of sessions (before averaging, if applicable).}
#'    \item{task_names}{The task part of the design matrix, after centering and scaling, but before any nuisance regression or prewhitening.}
#'  }
#'
# @importFrom ciftiTools read_cifti resample_gifti as.xifti remove_xifti
#' @import ciftiTools
#' @importFrom fMRItools unmask_mat dct_bases dct_convert match_input is_posNum
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
  dHRF=c(0, 1, 2),
  dHRF_as=c("auto", "nuisance", "task"),
  hpf=NULL,
  DCT=if(is.null(hpf)) {4} else {NULL},
  resamp_res=10000,
  # Below arguments shared with `BayesGLM`
  task_names = NULL,
  session_names = NULL,
  combine_sessions = TRUE,
  scale_BOLD = c("auto", "mean", "sd", "none"),
  scale_design = TRUE,
  Bayes = TRUE,
  #EM = FALSE,
  ar_order = 6,
  ar_smooth = 5,
  aic = FALSE,
  num.threads = 4,
  return_INLA = c("trimmed", "full", "minimal"),
  verbose = 1,
  meanTol = 1e-6,
  varTol = 1e-6#,emTol = 1e-3,
  ){

  EM <- FALSE
  emTol <- 1e-3

  # Preliminary steps. ---------------------------------------------------------
  ## Check simple arguments.
  ## These checks are in a separate function because they are shared with
  ## `BayesGLM_cifti`.
  argChecks <- BayesGLM_argChecks(
    combine_sessions = combine_sessions,
    scale_BOLD = scale_BOLD,
    scale_design = scale_design,
    Bayes = Bayes,
    EM = EM,
    ar_order = ar_order,
    ar_smooth = ar_smooth,
    aic = aic,
    num.threads = num.threads,
    return_INLA = return_INLA,
    verbose = verbose,
    meanTol = meanTol,
    varTol = varTol,
    emTol = emTol
  )
  scale_BOLD <- argChecks$scale_BOLD
  do_Bayesian <- argChecks$do_Bayesian
  do_EM <- argChecks$do_EM
  do_pw <- argChecks$do_pw
  return_INLA <- argChecks$return_INLA
  need_mesh <- do_Bayesian || (do_pw && ar_smooth > 0)

  # Brain structures.
  if ("both" %in% brainstructures) { brainstructures <- c("left", "right") }
  if ("all" %in% brainstructures) {
    message(
      "`brainstructures` is `all`, so using both left and right cortex. ",
      "Skipping subcortex (not implemented yet)."
    )
    brainstructures <- c("left","right") # "subcortical"
  }
  brainstructures <- fMRItools::match_input(
    brainstructures, c("left","right"),
    user_value_label="brainstructures"
  )
  do_left <- ('left' %in% brainstructures)
  do_right <- ('right' %in% brainstructures)

  # Nuisance arguments.
  dHRF <- as.numeric(match.arg(as.character(dHRF), c("0", "1", "2")))
  if (dHRF == 0) {
    if (identical(dHRF_as, "nuisance") || identical(dHRF_as, "task")) {
      warning("`dHRF_as` is only applicable if `dHRF > 0`. If `dHRF == 0`, there's no need to specify `dHRF_as`.")
    }
  }
  dHRF_as <- match.arg(dHRF_as, c("auto", "nuisance", "task"))
  if (!is.null(DCT)) {
    stopifnot(fMRItools::is_posNum(DCT, zero_ok=TRUE) && DCT==round(DCT))
    if (DCT==0) { DCT <- NULL }
  }
  if (!is.null(hpf)) {
    stopifnot(fMRItools::is_posNum(hpf, zero_ok=TRUE))
    if (hpf==0) { hpf <- NULL }
  }

  # xifti.
  #   Coerce to: a (length one) character vector, or a (length one) list of
  #   \code{"xifti"} objects.
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

  ## Sessions. -----------------------------------------------------------------
  # Name sessions and check compatibility of multi-session arguments
  nS <- nS_orig <- length(cifti_fname)
  if (nS==1) {
    combine_sessions <- FALSE
    if (is.null(session_names)) session_names <- 'single_session'
    if (!is.null(design)) design <- list(single_session = design)
    if (!is.null(onsets)) onsets <- list(single_session = onsets)
    if (!is.null(nuisance)) nuisance <- list(single_session = nuisance)
  } else {
    if (is.null(session_names)) session_names <- paste0('session', 1:nS)
    # if (length(session_names) == 1) { session_names <- paste0(session_names, 1:nsess) } # allow prefix?
    if (!is.null(design) && (length(design) != nS)) {
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`design` must be a list of length equal to the number of sessions ",
        " (or `NULL`, if onsets provided)."
      ))
    }
    if (!is.null(onsets) && (length(onsets) != nS)) {
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`onsets` must be a list of length equal to the number of sessions ",
        " (or `NULL`, if `design` provided)."
      ))
    }
    if (!is.null(nuisance) && (length(nuisance) != nS)) {
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`nuisance` must be a list of length equal to the number of sessions ",
        " (or `NULL`)."
      ))
    }
  }
  if (length(session_names) != nS) {
    stop('The length of `session_names` must match the number of sessions in `cifti_fname`.')
  }
  if(is.null(nuisance)) nuisance <- vector("list",length = nS)

  ## Surfaces: check or get. ---------------------------------------------------
  surf_list <- list(left=NULL, right=NULL)
  if (need_mesh) {
    if (do_left) {
      if (is.null(surfL_fname)) {
        if (verbose>0) cat("Using `ciftiTools` default inflated surface for the left cortex.\n")
        surfL_fname <- ciftiTools.files()$surf["left"]
      }
      surf_list$left <- read_surf(surfL_fname, resamp_res=resamp_res)
    }
    if (do_right) {
      if (is.null(surfR_fname)) {
        if (verbose>0) cat("Using `ciftiTools` default inflated surface for the right cortex.\n")
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

  ## `design`, `onsets`, `task_names`. -----------------------------------------
  ## Also, determine if we are doing multi-session modeling.
  if (!xor(is.null(design), is.null(onsets))) { stop('`design` or `onsets` must be provided, but not both.') }
  if (!is.null(design)) {
    do_multisesh <- inherits(design, "list")
    if (!do_multisesh) { stopifnot(inherits(design, "matrix") || inherits(design, "data.frame")) }
    d1 <- if (do_multisesh) { design[[1]] } else { design }
    task_names <- if (!is.null(task_names)) {
      task_names
    } else if (!is.null(names(d1))) {
      names(d1)
    } else {
      paste0("beta", seq(ncol(d1)))
    }
    rm(d1)
  }
  if (!is.null(onsets)) {
    if (is.null(TR)) { stop('Please provide `TR` if onsets provided') }
    do_multisesh <- inherits(onsets[[1]], "list")
    o1 <- if (do_multisesh) { onsets[[1]] } else { onsets }
    if (!do_multisesh) { stopifnot(inherits(o1, "matrix") || inherits(o1, "data.frame")) }
    task_names <- if (!is.null(task_names)) {
      task_names
    } else if (!is.null(names(o1))) {
      names(o1)
    } else {
      paste0("beta", seq(length(o1)))
    }
    rm(o1)
  }

  # Data setup. ----------------------------------------------------------------
  if (verbose>0) cat('Setting up data:\n')

  ## xifti things. -------------------------------------------------------------
  ### For each session, separate the CIFTI data into left/right/sub and read in files
  BOLD_list <- list(left=NULL, right=NULL)
  mwallL <- mwallR <- NULL
  ntime <- vector("numeric", nS)

  for (ss in seq(nS)) {
    if (nS > 1) if (verbose>0) cat(paste0('\tReading in data for session ', ss,'.\n'))

    if (is_xifti) {
      xii_ss <- cifti_fname[[ss]]
      xii_ss_res <- infer_resolution(xii_ss)
      if (!is.null(resamp_res) && any(infer_resolution(xii_ss)!=resamp_res)) {
        xii_ss <- resample_xifti(xii_ss, resamp_res=resamp_res)
      }
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
      if (is.null(xii_ss$data$cortex_left)) {
        stop("No left cortex data for session ", ss, ".")
      }
      BOLD_list[["left"]][[ss]] <- fMRItools::unmask_mat(xii_ss$data$cortex_left, mwallL)
    }
    if (do_right) {
      if (is.null(xii_ss$data$cortex_right)) {
        stop("No right cortex data for session ", ss, ".")
      }
      BOLD_list[["right"]][[ss]] <- fMRItools::unmask_mat(xii_ss$data$cortex_right, mwallR)
    }
  }

  BOLD_list <- BOLD_list[!vapply(BOLD_list, is.null, FALSE)]
  if (need_mesh) { surf_list <- surf_list[!vapply(surf_list, is.null, FALSE)] }

  ## Design and nuisance matrices. ---------------------------------------------
  if (is.null(design)) {
    if (verbose>0) cat("\tMaking design matrices.\n")
    design <- vector("list", nS)

    for (ss in seq(nS)) {
      HRF_ss <- make_HRFs(
        onsets[[ss]], TR=TR, duration=ntime[ss],
        dHRF=dHRF, dHRF_as=dHRF_as,
        verbose=ss==1
      )
      design[[ss]] <- HRF_ss$design
      if (!is.null(HRF_ss$nuisance)) {
        if (!is.null(nuisance[[ss]])) {
          nuisance[[ss]] <- cbind(nuisance[[ss]], HRF_ss$nuisance)
        } else {
          nuisance[[ss]] <- HRF_ss$nuisance
        }
      }
    }
  }

  # Check that design matrix names are consistent across sessions.
  if (nS > 1) {
    tmp <- sapply(design, colnames)
    if(length(task_names) == 1) {
      num_names <- length(unique(tmp))
      if (num_names > 1) stop('task names must match across sessions for multi-session modeling')
    } else {
      num_names <- apply(tmp, 1, function(x) length(unique(x))) #number of unique names per row
      if (max(num_names) > 1) stop('task names must match across sessions for multi-session modeling')
    }
  }

  # Warn the user if the number of design matrix columns exceeds five.
  if (Bayes && ncol(design[[1]]) > 5) {
    message("The number of design matrix columns exceeds five. INLA computation may be very slow. To avoid stalling, you can quit this function call now and modify the analysis. For example, model some signals as nuisance rather than tasks.")
    Sys.sleep(10)
  }

  task_names <- colnames(design[[1]]) # because if dHRF > 0, there will be more task_names.

  # Scale design matrix. (Here, rather than in `BayesGLM`, b/c it's returned.)
  design <- if (scale_design) {
    sapply(design, scale_design_mat, simplify = F)
  } else {
    sapply(design, scale, scale = F, simplify = F)
  }

  # Add DCT bases.
  DCTs <- vector("numeric", nS)
  for (ss in 1:nS) {
    # DCT highpass filter
    if (!is.null(hpf) || !is.null(DCT)) {
      # Get the num. of bases for this session.
      if (!is.null(hpf)) {
        DCTs[ss] <- round(fMRItools::dct_convert(ntime[ss], TR, f=hpf))
      } else {
        DCTs[ss] <- DCT
      }
      # Generate the bases and add them.
      DCTb_ss <- fMRItools::dct_bases(ntime[ss], DCTs[ss])
      if (DCTs[ss] > 0) {
        if (!is.null(nuisance[[ss]])) {
          nuisance[[ss]] <- cbind(nuisance[[ss]], DCTb_ss)
        } else {
          nuisance[[ss]] <- DCTb_ss
        }
      }
    }
  }

  # Do GLM. --------------------------------------------------------------------
  BayesGLM_results <- list(left = NULL, right = NULL)

  # >> Loop through brainstructures to complete the analyses on the different hemispheres ----
  for (bb in brainstructures) {

    if (verbose>0) cat(paste0(toupper(bb)," cortex analysis:\n"))

    # set up session list
    session_data <- vector('list', nS)
    names(session_data) <- session_names
    for (ss in seq(nS)) {
      session_data[[ss]] <- list(
        BOLD = t(BOLD_list[[bb]][[ss]]),
        design = design[[ss]]
      )
      if (!is.null(nuisance[[ss]])) {
        session_data[[ss]]$nuisance <- nuisance[[ss]]
      }
    }

    BayesGLM_out <- BayesGLM(
      data = session_data,
      vertices = surf_list[[bb]]$vertices,
      faces = surf_list[[bb]]$faces,
      mesh = NULL,
      mask = NULL,
      task_names = NULL, # in `session_data`
      session_names = session_names,
      combine_sessions = combine_sessions,
      scale_BOLD = scale_BOLD,
      scale_design = FALSE, # done above
      Bayes = do_Bayesian,
      #EM = do_EM,
      ar_order = ar_order,
      ar_smooth = ar_smooth,
      aic = aic,
      num.threads = num.threads,
      return_INLA = return_INLA,
      verbose = verbose,
      meanTol=meanTol,
      varTol=varTol#,emTol=emTol,
    )

    BayesGLM_results[[bb]] <- BayesGLM_out

    # update session info if averaged over sessions
    if (bb == brainstructures[1] && combine_sessions) {
      session_names <- BayesGLM_out$session_names
      nS <- 1
    }

    rm(BayesGLM_out); gc()
  }
  names(BayesGLM_results)[names(BayesGLM_results)=="left"] <- "cortex_left"
  names(BayesGLM_results)[names(BayesGLM_results)=="right"] <- "cortex_right"

  ### CONSTRUCT BETA ESTIMATES AS CIFTI OBJECTS
  if (verbose>0) cat("Formatting results.\n")
  task_cifti_classical <- task_cifti <- vector('list', nS)
  names(task_cifti_classical) <- names(task_cifti) <- session_names
  datL <- datR <- NULL
  for (ss in seq(nS)) {

    # CLASSICAL GLM
    if (do_left) {
      datL <- BayesGLM_results$cortex_left$result_classical[[ss]]$estimates
      mwallL <- !is.na(datL[,1]) # update b/c mask2 can change the medial wall
      datL <- datL[mwallL,]
    }
    if (do_right) {
      datR <- BayesGLM_results$cortex_right$result_classical[[ss]]$estimates
      mwallR <- !is.na(datR[,1])
      datR <- datR[mwallR,]
    }
    task_cifti_classical[[ss]] <- as.xifti(
      cortexL = datL,
      cortexL_mwall = mwallL,
      cortexR = datR,
      cortexR_mwall = mwallR
    )
    task_cifti_classical[[ss]]$meta$cifti$names <- task_names

    # BAYESIAN GLM
    if (do_Bayesian) {
      if (do_left) {
        datL <- BayesGLM_results$cortex_left$task_estimates[[ss]]
        mwallL <- !is.na(datL[,1])
        datL <- datL[mwallL,]
      }
      if (do_right) {
        datR <- BayesGLM_results$cortex_right$task_estimates[[ss]]
        mwallR <- !is.na(datR[,1])
        datR <- datR[mwallR,]
      }
      task_cifti[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = mwallL,
        cortexR = datR,
        cortexR_mwall = mwallR
      )
      task_cifti[[ss]]$meta$cifti$names <- task_names
    }
  }

  result <- list(
    task_estimates_xii = list(
      Bayes = task_cifti,
      classical = task_cifti_classical
    ),
    task_names = task_names,
    session_names = session_names,
    n_sess_orig = nS_orig,
    # task part of design matrix after centering/scaling but
    #   before nuisance regression and prewhitening.
    design = design,
    BayesGLM_results = BayesGLM_results
  )

  if (verbose>0) cat('Done!\n')

  class(result) <- "BayesGLM_cifti"
  result
}
