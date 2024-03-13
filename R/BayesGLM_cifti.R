#' BayesGLM for CIFTI
#'
#' Performs spatial Bayesian GLM on the cortical surface for task fMRI
#'  activation.
#'
#' @section INLA latent fields limit:
#'  INLA computation times increase greatly when the number of columns in the
#'  design matrix exceeds five. So if there are more than five tasks, or
#'  three or more tasks each with its temporal derivative being modeled as a
#'  task, \code{BayesGLM} will raise a warning. In cases like the latter, we
#'  recommend modeling the temporal derivatives as nuisance signals using the
#'  \code{nuisance} argument, rather than modeling them as fields.
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
#'
#'  If \code{cifti_fname} is a \code{"xifti"} object or list of \code{"xifti"}
#'  objects, its surfaces, if any, will be used for Bayesian modeling. However,
#'  \code{surfL_fname} and \code{surfR_fname}, if provided, will override any
#'  included surfaces.
#' @param surfL_fname,surfR_fname Left or right cortex surface geometry in
#'  GIFTI format ("*.surf.gii"). This can be a file path to a GIFTI file or a
#'  \code{"surf"} object from the \code{ciftiTools} package. This argument is
#'  only used if \code{brainstructures} includes the corresponding hemisphere
#'  and \code{Bayes==TRUE}. If it's not provided, and there are no surfaces in
#'  \code{cifti_fname}, the HCP group-average inflated surface included in the
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
#'   represents the expected BOLD response for each field, a convolution of
#'   the hemodynamic response function (HRF) and the task stimulus. Field names
#'   should be the column names. For multi-session modeling, this argument should be a list of
#'   such matrices, and column names must match and align across sessions. If any field is
#'   missing from any session, include a column of zeros in its place. HRF derivatives
#'   may be included in \code{design} to model them spatially, or in \code{nuisance}
#'   to treat them as nuisance signals. Note that INLA computation times increase
#'   if the design matrix has more than five columns.
#'
#'   \code{onsets} is a list where each element represents a task, and the value of each element is a
#'   matrix of onsets (first column) and durations (second column) for each
#'   stimuli (each row) of the corresponding task. The units of both columns
#'   is seconds. For multi-session modeling, this argument should be a list of
#'   such lists, and the task names must match and align across sessions. If any task is
#'   missing from any session, include a list element equal to NA in its place.
#'
#'   \code{TR} is the temporal resolution of the data, in seconds. Required if \code{onsets} provided.
#'
#' @param design_multiple (Optional) A \eqn{T \times K \times D} array of \eqn{D}
#' different design matrices for model comparison.  If provided, onsets and design will be ignored.
#' TO DO: Allow differing numbers of regressors across D models, pad with NAs or 0.
# @param task_names Names of tasks represented in design matrix.  For multi-session
#' modeling, this must be provided if not all tasks are present for all sessions.
#' @inheritParams session_names_Param
#' @param nuisance (Optional) A \eqn{T \times J} matrix of nuisance signals.
#'  These are regressed from the fMRI data and the design matrix prior to the
#'  GLM computation. For multi-session modeling, this argument should be a list
#'  of such matrices.
#' @param dHRF,dHRF_as Only applicable if \code{onsets} and \code{TR} are
#'  provided. These arguments enable the modeling of HRF derivatives.
#'
#'  Set \code{dHRF} to \code{1} to model the temporal derivative of the HRF
#'  (default), \code{2} to add the dispersion derivative too, or \code{0} to
#'  include only the main HRF regressor.
#'
#'  If \code{dHRF > 0}, \code{dHRF_as} controls whether the derivatives are
#'  modeled as \code{"nuisance"} signals to regress out, \code{"field"}, or
#'  \code{"auto"} (default) to treat them as fields unless the total
#'  number of columns in the design matrix would exceed five (for computational
#'  reasons).
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
#' @param nbhd_order For volumetric data, what order neighborhood around data
#' locations to keep? (0 = no neighbors, 1 = 1st-order neighbors, 2 = 1st- and
#' 2nd-order neighbors, etc.). Smaller values will provide greater computational
#' efficiency at the cost of higher variance around the edge of the data.
#' @param buffer For volumetric data, size of extra voxels layers around the
#' bounding box, in terms of voxels. Set to NULL for no buffer. (Recommended not
#' to change unless you know what you're doing. Instead to reduce the number of
#' boundary voxels, adjust \code{nbhd_order}).
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
# @inheritParams combine_sessions_Param
#' @param meanTol,varTol Tolerance for mean and variance of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Default: \code{1e-6} for both.
# @inheritParams emTol_Param
#' @inheritParams trim_INLA_Param
#'
#' @return An object of class \code{"BayesGLM_cifti"}: a list with elements
#'  \describe{
#'    \item{betas_Bayesian}{The field coefficients for the Bayesian model.}
#'    \item{betas_classical}{The field coefficients for the classical model.}
#'    \item{GLMs_Bayesian}{The entire list of GLM results, except for parameters estimated for the classical model.}
#'    \item{GLMs_classical}{Parameters estimated for the classical model from the GLM.}
#'    \item{session_names}{The names of the sessions.}
#'    \item{n_sess_orig}{The number of sessions (before averaging, if applicable).}
#'    \item{field_names}{Column names of the fields in the design matrix.}
#'  }
#'
# @importFrom ciftiTools read_cifti resample_gifti as.xifti remove_xifti
#' @import ciftiTools
#' @importFrom fMRItools unmask_mat dct_bases dct_convert match_input is_posNum
#' @importFrom matrixStats rowVars rowSums2 colVars
#' @importFrom parallel detectCores
#' @importFrom Matrix bdiag
#'
#' @export
BayesGLM_cifti <- function(
  cifti_fname,
  surfL_fname=NULL,
  surfR_fname=NULL,
  brainstructures=c('left','right'),
  design=NULL,
  design_multiple=NULL,
  onsets=NULL,
  TR=NULL,
  #task_names = NULL, #disabled this, user must provide through design or onsets
  session_names = NULL,
  nuisance=NULL,
  dHRF=c(1, 0, 2),
  dHRF_as=c("auto", "nuisance", "field"),
  hpf=NULL,
  DCT=if(is.null(hpf)) {4} else {NULL},
  resamp_res=10000,
  nbhd_order=1, buffer=c(1,1,3,4,4),
  # Below arguments shared with `BayesGLM`.
  #combine_sessions = TRUE,
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
    #combine_sessions = combine_sessions,
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

  # Brain structures.
  if ("both" %in% brainstructures) { brainstructures <- c("left", "right") }
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }
  brainstructures <- match_input(
    brainstructures, c("left","right","subcortical"),
    user_value_label="brainstructures"
  )
  do_left <- ('left' %in% brainstructures)
  do_right <- ('right' %in% brainstructures)
  do_sub <- ('subcortical' %in% brainstructures)
  do_cort <- do_left || do_right

  #need_mesh <- do_Bayesian #always need meshes for Bayesian modeling
  #if(do_pw && ar_smooth > 0 && do_cort) need_mesh <- TRUE #also need meshes for ar-smoothing if doing cortical modeling

  # # Temporary
  # if (need_mesh && do_sub) {
  #   stop("Bayesian modeling and AR smoothing both require spatial modeling, which is currently not availble for the subcortex.")
  # }

  # Nuisance arguments.
  dHRF <- as.numeric(match.arg(as.character(dHRF), c("1", "0", "2")))
  if (dHRF == 0) {
    if (identical(dHRF_as, "nuisance") || identical(dHRF_as, "field")) {
      warning("`dHRF_as` is only applicable if `dHRF > 0`. If `dHRF == 0`, there's no need to specify `dHRF_as`.")
    }
  }
  dHRF_as <- match.arg(dHRF_as, c("auto", "nuisance", "field"))
  if (!is.null(DCT)) {
    stopifnot(is_posNum(DCT, zero_ok=TRUE) && DCT==round(DCT))
    if (DCT==0) { DCT <- NULL }
  }
  if (!is.null(hpf)) {
    stopifnot(is_posNum(hpf, zero_ok=TRUE))
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

  # Use the `"xifti"` surfaces, if available.
  if (is_xifti) {
    if (is.null(surfL_fname) && !is.null(cifti_fname[[1]]$surf$cortex_left)) {
      all_left_surfs <- c(
        lapply(cifti_fname, function(q){q$surf$cortex_left}), list(NULL)
      )
      if (length(unique(all_left_surfs)) > 2) {
        warning("Using left surface from the first `xifti` for all modeling. Ignoring the other left surfaces.\n")
      } else {
        cat("Using left surface from the `xifti` data.\n")
      }
      surfL_fname <- cifti_fname[[1]]$surf$cortex_left
    }
    if (is.null(surfR_fname) && !is.null(cifti_fname[[1]]$surf$cortex_right)) {
      all_right_surfs <- c(
        lapply(cifti_fname, function(q){q$surf$cortex_right}), list(NULL)
      )
      if (length(unique(all_right_surfs)) > 2) {
        warning("Using right surface from the first `xifti` for all modeling. Ignoring the other right surfaces.\n")
      } else {
        cat("Using right surface from the `xifti` data.\n")
      }
      surfR_fname <- cifti_fname[[1]]$surf$cortex_right
    }
  }

  # [TO DO]: If input is a `"xifti"`, infer `resamp_res`
  # or maybe just add surfaces to the `"xifti"` using `add_surf` and that will handle the
  # difference in resolution.

  ## Sessions. -----------------------------------------------------------------
  # Name sessions and check compatibility of multi-session arguments
  nS <- nS_orig <- length(cifti_fname)

  if (!is.null(session_names) && (length(session_names) != nS))
    stop('The length of `session_names` must match the number of sessions in `cifti_fname`.')

  if (nS==1) {
    if (verbose>0) { cat("Preparing to analyze a single task fMRI session.\n") }

    #combine_sessions <- FALSE
    if (is.null(session_names)) session_names <- 'single_session'
    if (!is.null(design)) design <- list(single_session = design)
    if (!is.null(design_multiple)) design_multiple <- list(single_session = design_multiple)
    if (!is.null(onsets)) onsets <- list(single_session = onsets)
    if (!is.null(nuisance)) nuisance <- list(single_session = nuisance)

  } else {

    if (verbose>0) { cat(paste0(
      "Preparing to analyze ",nS,
      " task fMRI sessions with a common set of tasks.\n"
    )) }

    #name sessions
    if (is.null(session_names)) session_names <- paste0('session', 1:nS)

    #check that length of design = nS
    if (!is.null(design) && (length(design) != nS)){
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`design` must be a list of length equal to the number of sessions ",
        " (or `NULL`, if `onsets` or `design_multiple` provided)."
      ))
    }

    #check that length of design_multiple = nS
    if (!is.null(design_multiple) && (length(design_multiple) != nS)){
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`design_multiple` must be a list of length equal to the number of sessions ",
        " (or `NULL`, if `onsets` or `design` provided)."
      ))
    }

    #check that length of onsets = nS
    if (!is.null(onsets) && (length(onsets) != nS)) {
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`onsets` must be a list of length equal to the number of sessions ",
        " (or `NULL`, if `design` or `design_multiple` provided)."
      ))
    }

    #check that length of nuisance = nS
    if (!is.null(nuisance) && (length(nuisance) != nS)) {
      stop(paste(
        "If multiple sessions provided (because `cifti_fname` is a vector), ",
        "`nuisance` must be a list of length equal to the number of sessions ",
        " (or `NULL`)."
      ))
    }
  }

  if (!is.null(nuisance)) {
    stopifnot(all(vapply(nuisance, is.matrix, FALSE)))
    stopifnot(all(vapply(nuisance, is.numeric, FALSE)))
  } else {
    nuisance <- vector("list", length = nS)
  }

  ## Surfaces/SPDE: check or get. ---------------------------------------------------
  spatial_list <- list(left=NULL, right=NULL, subcortical=NULL)
  if (do_left) {
    if (is.null(surfL_fname)) {
      if (verbose>0) cat("Using `ciftiTools` default inflated fs_LR surface for the left cortex.\n")
      surfL_fname <- ciftiTools.files()$surf["left"]
    }
    if (suppressMessages(is.surf(surfL_fname))) {
      if (!is.null(resamp_res)) { surfL_fname <- resample_surf(surfL_fname, resamp_res=resamp_res)  }
    } else {
      surfL_fname <- read_surf(surfL_fname, resamp_res=resamp_res)
    }
    spatial_list$left <- surfL_fname
  }
  if (do_right) {
    if (is.null(surfR_fname)) {
      if (verbose>0) cat("Using `ciftiTools` default inflated fs_LR surface for the right cortex.\n")
      surfR_fname <- ciftiTools.files()$surf["right"]
    }
    if (suppressMessages(is.surf(surfR_fname))) {
      if (!is.null(resamp_res)) { surfR_fname <- resample_surf(surfR_fname, resamp_res=resamp_res)  }
    } else {
      surfR_fname <- read_surf(surfR_fname, resamp_res=resamp_res)
    }
    spatial_list$right <- surfR_fname  }
  if(do_sub) {
    #get mask and labels from xii data
    if (is_xifti) xii_1 <- cifti_fname[[1]] else xii_1 <- read_cifti(cifti_fname[1], brainstructures='sub', idx=1)
    mask <- xii_1$meta$subcort$mask
    labels <- xii_1$meta$subcort$labels
    labels_img <- mask*0; labels_img[mask==TRUE] <- labels
    spatial_list$subcortical <- labels_img #this is a 3D array of labels
  }

  ## `design`, `onsets`, `task_names`. -----------------------------------------

  if(!is.null(design)) { task_names <- colnames(design[[1]]); nK <- ncol(design[[1]]) } #task_names could still be NULL
  if(!is.null(onsets)) { task_names <- names(onsets[[1]]); nK <- length(onsets[[1]]) } #task_names could still be NULL
  if(!is.null(design_multiple)) { task_names <- field_names <- colnames(design_multiple[[1]][,,1]); nK <- dim(design_multiple[[1]])[2] } #task_names could still be NULL
  if(is.null(task_names)) stop('Task/field names must be specified through `onsets` or `design`. See documentation for details.')


  if(is.null(design_multiple)){

    if (!xor(is.null(design), is.null(onsets))) { stop('`design` or `onsets` must be provided, but not both.') }

    ### Case 1: Design matrix provided
    if (!is.null(design)) {
      #check format of design
      if(!all(vapply(design, function(q){is.matrix(q) | is.data.frame(q)}, FALSE))) {
        stop('`design` must be a numeric TxK matrix, or list of such matrices for multi-session analysis.')
      }
      if(!all(vapply(design, is.numeric, FALSE))) {
        stop('`design` must be a numeric TxK matrix, or list of such matrices for multi-session analysis.')
      }
      if(!all(vapply(design, function(q){ncol(q)>0}, FALSE))) {
        stop('`design` has data with zero columns. Please fix.')
      }
    }

    ### Case 2: Onsets and TR provided
    for (ss in seq(nS)) {
      onsets_ss <- onsets[[ss]]
      if (!is.null(onsets_ss)) {
        if (is.null(TR)) { stop('Please provide `TR` if `onsets` provided.') }
        #check format of onsets_ss
        if(!all(vapply(onsets_ss, is_onsets, FALSE))) {
          if(nS > 1) cat(paste0('For session number ',ss,', '))
          stop('`onsets` is not a valid format.  Each list element must be a non-empty numeric matrix or data frame.')
        }
      }
    }

  } else {

    ### Case 3: Multiple design matrices provided

    #check format of design_multiple
    dclass <- sapply(design_multiple, function(x){
      if(inherits(x, "array") && length(dim(x))==3) return("array") else return(NA)
    })
    if(any(is.na(dclass))) stop('`design_multiple` must be a 3D array, or list of 3D arrays for multi-session analysis')
  }



  # Data setup. ----------------------------------------------------------------
  if (verbose>0) cat('Setting up data:\n')

  ## xifti things. -------------------------------------------------------------
  ### For each session, separate the CIFTI data into left/right/sub and read in files
  BOLD_list <- list(left=NULL, right=NULL, subcortical=NULL)
  mwallL <- mwallR <- NULL
  ntime <- vector("numeric", nS)

  for (ss in seq(nS)) {
    if (verbose>0) {
      if (nS==1) {
        cat('\tReading BOLD data.\n')
      } else {
        if (ss==1) { cat('\tReading BOLD data for session ') }
        if (ss!=nS) { cat(paste0(ss, ", ")) } else { cat(paste0(ss, ".\n")) }
      }
    }

    if (is_xifti) {
      xii_ss <- cifti_fname[[ss]]
      xii_ss_res <- ciftiTools::infer_resolution(xii_ss)
      if(do_left | do_right){ #only do resampling for surface data
        if (!is.null(resamp_res) && any(ciftiTools::infer_resolution(xii_ss)!=resamp_res)) {
          xii_ss <- resample_xifti(xii_ss, resamp_res=resamp_res)
        }
      }
    } else {
      xii_ss <- read_cifti(
        cifti_fname[ss], brainstructures=brainstructures,
        resamp_res=resamp_res
      )
    }

    mwallL_ss <- xii_ss$meta$cortex$medial_wall_mask$left
    mwallR_ss <- xii_ss$meta$cortex$medial_wall_mask$right
    submeta_ss <- if (do_sub) { xii_ss$meta$subcort } else { NULL }
    ntime[ss] <- ncol(xii_ss)

    # Get medial wall mask, or check that it matches.
    if (ss == 1) {
      if (do_left) { mwallL <- mwallL_ss }
      if (do_right) { mwallR <- mwallR_ss }
      # [TO DO] Check compatibility with `spatial_list`
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
      BOLD_list[["left"]][[ss]] <- unmask_mat(xii_ss$data$cortex_left, mwallL)
    }
    if (do_right) {
      if (is.null(xii_ss$data$cortex_right)) {
        stop("No right cortex data for session ", ss, ".")
      }
      BOLD_list[["right"]][[ss]] <- unmask_mat(xii_ss$data$cortex_right, mwallR)
    }
    if (do_sub) {
      if (is.null(xii_ss$data$subcort)) {
        stop("No subcortical data for session ", ss, ".")
      }
      BOLD_list[["subcortical"]][[ss]] <- xii_ss$data$subcort
    }
  }

  #remove un-needed brain structures
  BOLD_list <- BOLD_list[!vapply(BOLD_list, is.null, FALSE)]
  spatial_list <- spatial_list[!vapply(spatial_list, is.null, FALSE)]

  ## Design and nuisance matrices. ---------------------------------------------

  if (!is.null(onsets)) {

    if (verbose>0) cat("\tConstructing design matrix based on onsets provided.\n")

    #determine whether to model HRF derivatives as field or nuisance if "auto"
    nK <- length(onsets[[1]])
    if (dHRF > 0 && dHRF_as=="auto") {
      nJ <- (dHRF+1) * nK # number of design matrix columns
      if (nJ > 5) {
        if (verbose) {
          message("Treating the HRF derivatives as nuisance signals.")
        }
        dHRF_as <- "nuisance"
      } else {
        if (verbose) {
          message("Modeling the HRF derivatives as fields.")
        }
        dHRF_as <- "field"
      }
    }

    #initialize lists
    design <- design_FIR <- vector("list", nS)
    stimulus <- HRFs <- FIR <- vector("list", nS)

    for (ss in seq(nS)) {

      ### Construct stimulus function, HRFs, derivatives, and FIR basis set

      HRF_ss <- make_HRFs(
        onsets[[ss]],
        TR=TR,
        nT=ntime[ss],
        dHRF=2 # Damon: default `dHRF` value for `make_HRFs` used to be `2`. And we were using this default value. Check if `2` is necessary here still, or if it can be what `dHRF` is for `BayesGLM_cifti`.
      )
      stimulus[[ss]] <- HRF_ss$stimulus #TxK matrix
      HRFs[[ss]] <- HRF_ss$HRF_convolved  #TxKx3 array -- allocate 2nd and 3rd dims to design or nuisance
      FIR[[ss]] <- HRF_ss$FIR #T x K x nFIR -- just return this for now, don't use in modeling
      task_names_ss <- colnames(HRF_ss$stimulus)
      if(ss==1) HRF <- HRF_ss$HRF

      if(any(task_names_ss != task_names)) stop(paste0('Tasks in session ',nS,' in `onsets` do not match other sessions. Please input NA for any missing tasks, as described in function documentation.'))

      ### Construct design and nuisance matrices from main HRF and derivatives

      design[[ss]] <- HRF_ss$design #main HRF
      if(dHRF > 0){
        for(dd in 1:dHRF){
          dname_dd <- switch(dd, "dHRF", "ddHRF")
          dHRF_sd <- as.matrix(HRF_ss$HRF_convolved[,,(dd+1)], ncol=nK)
          colnames(dHRF_sd) <- paste0(task_names_ss, "_", dname_dd) # now these are field names
          if(dHRF_as == "nuisance") nuisance[[ss]] <- cbind(nuisance[[ss]], dHRF_sd) #if nuisance[[ss]] is NULL this still works
          if(dHRF_as == "field") design[[ss]] <- cbind(design[[ss]], dHRF_sd)
        }
      }

      ### Construct FIR design matrix
      if(!is.null(HRF_ss$FIR)){
        for(kk in 1:nK){
          FIR_kk <- HRF_ss$FIR[,kk,]
          colnames(FIR_kk) <- paste0(task_names[kk], '_FIR',1:ncol(FIR_kk))
          design_FIR[[ss]] <- cbind(design_FIR[[ss]], FIR_kk)
        }
      }

    } #end loop over sessions
  } #end design matrix construction from onsets

  if(!is.null(design)){

    # Check that design matrix names are consistent across sessions.
    if (nS > 1) {
      tmp <- lapply(design, colnames) #list of column names for all sessions
      tmp2 <- sapply(tmp, function(x) all.equal(x, tmp[[1]])) #vector with TRUE for sessions whose names match session 1
      if(!all(tmp2 == TRUE)) stop(paste0('Fields do not match across sessions in `design`. Please input empty values for any missing tasks, as described in function documentation.'))
    }

    # Warn the user if the number of design matrix columns exceeds five.
    if (Bayes && ncol(design[[1]]) > 5) {
      message("The number of regressors to be modeled spatially exceeds five. INLA computation may be slow. Consider reducing the number of design matrix columns, e.g. by modeling HRF derivatives as nuisance")
      Sys.sleep(10)
    }

    field_names <- colnames(design[[1]]) # because if dHRF > 0, there will be more design matrix columns

    # Scale design matrix. (Here, rather than in `BayesGLM`, b/c it's returned.)
    design <- if(scale_design) {
      sapply(design, scale_design_mat, simplify = F)
    } else {
      sapply(design, scale, scale = F, simplify = F)
    }

    stimulus <- HRFs <- FIR <- design_FIR <- vector("list", nS) #leave empty in this case

  }

  # [TO DO]
  #field_names <- field_names

  # Add DCT bases to nuisance matrix

  DCTs <- vector("numeric", nS)
  for (ss in 1:nS) {
    # DCT highpass filter
    if (!is.null(hpf) || !is.null(DCT)) {
      # Get the num. of bases for this session.
      if (!is.null(hpf)) {
        DCTs[ss] <- round(dct_convert(ntime[ss], TR, f=hpf))
        if(verbose > 0 & ss==1) cat(paste('\tIncluding',DCTs[ss],'DCT basis functions in the model based for hpf =',hpf,'Hz\n'))
      } else {
        DCTs[ss] <- DCT
      }
      # Generate the bases and add them.
      DCTb_ss <- dct_bases(ntime[ss], DCTs[ss])
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
  BayesGLM_results <- list(left = NULL, right = NULL, subcortical = NULL)

  if(!is.null(design_multiple)){
    design <- lapply(design_multiple, function(x) x[,,1]) #use first model from each session as a placeholder in session_data below
  }

  # >> Loop through brainstructures ----
  for (bb in brainstructures) {
    bname <- list(left="left cortex", right="right cortex", subcortical="subcortical")[[bb]]

    if (verbose>0) cat(paste0(toupper(bname)," analysis:\n"))

    # set up session list
    session_data <- vector('list', nS)
    names(session_data) <- session_names
    for (ss in seq(nS)) {
      session_data[[ss]] <- list(
        BOLD = t(BOLD_list[[bb]][[ss]]),
        design = design[[ss]]
        #design_FIR = design_FIR[[ss]]
      )
      if (!is.null(nuisance[[ss]])) {
        session_data[[ss]]$nuisance <- nuisance[[ss]]
      }
    }

    #Extract spatial information and pass to BayesGLM
    vertices <- faces <- labels <- NULL
    if(do_cort){
      vertices <- spatial_list[[bb]]$vertices
      faces <- spatial_list[[bb]]$faces
    } else {
      labels <- spatial_list[[bb]]
    }

    #Set scale_design = FALSE if not fitting multiple models
    if(is.null(design_multiple)) scale_design <- FALSE # scaling already done above

    BayesGLM_out <- BayesGLM(
      data = session_data,
      design_multiple = design_multiple,
      vertices = vertices,
      faces = faces,
      labels = labels,
      nbhd_order = nbhd_order,
      buffer = buffer,
      scale_design = scale_design,
      Bayes = do_Bayesian,
      #combine_sessions = combine_sessions,
      scale_BOLD = scale_BOLD,
      #EM = do_EM,
      ar_order = ar_order,
      ar_smooth = ar_smooth,
      aic = aic,
      num.threads = num.threads,
      return_INLA = return_INLA,
      verbose = verbose,
      meanTol = meanTol,
      varTol = varTol#,emTol=emTol,
    )

    BayesGLM_results[[bb]] <- BayesGLM_out

    # # update session info if averaged over sessions
    # if (bb == brainstructures[1] && combine_sessions) {
    #   session_names <- BayesGLM_out$session_names
    #   nS <- 1
    # }

    rm(BayesGLM_out); gc()
  }
  names(BayesGLM_results)[names(BayesGLM_results)=="left"] <- "cortex_left"
  names(BayesGLM_results)[names(BayesGLM_results)=="right"] <- "cortex_right"

  ### CONSTRUCT BETA ESTIMATES AS CIFTI OBJECTS

  if (verbose>0) cat("Formatting results.\n")

  field_cifti_classical <- field_cifti <- vector('list', nS)
  names(field_cifti_classical) <- names(field_cifti) <- session_names

  if(is.null(design_multiple)){

    datL <- datR <- datSub <- NULL
    for (ss in seq(nS)) {

      # CLASSICAL GLM
      if (do_left) {
        datL <- BayesGLM_results$cortex_left$result_classical[[ss]]$estimates
        mwallL <- !is.na(datL[,1]) # update b/c mask2 can change the medial wall
        datL <- datL[mwallL,]
        colnames(datL) <- NULL
      }
      if (do_right) {
        datR <- BayesGLM_results$cortex_right$result_classical[[ss]]$estimates
        mwallR <- !is.na(datR[,1])
        datR <- datR[mwallR,]
        colnames(datR) <- NULL
      }
      if (do_sub) {
        datSub <- BayesGLM_results$subcort$result_classical[[ss]]$estimates
        colnames(datSub) <- NULL
      }
      field_cifti_classical[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = mwallL,
        cortexR = datR,
        cortexR_mwall = mwallR,
        subcortVol = datSub,
        subcortLabs = submeta_ss$labels,
        subcortMask = submeta_ss$mask
      )
      field_cifti_classical[[ss]]$meta$subcort$trans_mat <- submeta_ss$trans_mat
      field_cifti_classical[[ss]]$meta$subcort$trans_units <- submeta_ss$trans_units
      field_cifti_classical[[ss]]$meta$cifti$names <- field_names

      # BAYESIAN GLM
      if (do_Bayesian) {
        if (do_left) {
          datL <- BayesGLM_results$cortex_left$field_estimates[[ss]]
          mwallL <- !is.na(datL[,1])
          datL <- datL[mwallL,]
          colnames(datL) <- NULL
        }
        if (do_right) {
          datR <- BayesGLM_results$cortex_right$field_estimates[[ss]]
          mwallR <- !is.na(datR[,1])
          datR <- datR[mwallR,]
          colnames(datR) <- NULL
        }
        if (do_sub) {
          datSub <- BayesGLM_results$subcortical$field_estimates[[ss]]
          colnames(datSub) <- NULL
        }
        field_cifti[[ss]] <- as.xifti(
          cortexL = datL,
          cortexL_mwall = mwallL,
          cortexR = datR,
          cortexR_mwall = mwallR,
          subcortVol = datSub,
          subcortLabs = submeta_ss$labels,
          subcortMask = submeta_ss$mask
        )
        field_cifti[[ss]]$meta$subcort$trans_mat <- submeta_ss$trans_mat
        field_cifti[[ss]]$meta$subcort$trans_units <- submeta_ss$trans_units
        field_cifti[[ss]]$meta$cifti$names <- field_names
      }
    }

    names(design) <- names(nuisance) <- names(stimulus) <- session_names
    #if(length(FIR) > 0) names(FIR) <- names(design_FIR) <- session_names #this check doesn't work for multi-session modeling
    result_multiple_xii <- NULL

  } else {

    #format as xii: (1) index of best model and (2) locations of no signal
    result_multiple_xii <- vector('list', length=2)
    names(result_multiple_xii) <- c('bestmodel','pval_F')
    for(meas in c('bestmodel','pval_F')){

      datL <- datR <- datSub <- NULL

      if (do_left) {
        datL <- BayesGLM_results$cortex_left$result_multiple[[meas]]*1
        mwallL <- !is.na(datL) # update b/c mask2 can change the medial wall
        datL <- datL[mwallL]
      }
      if (do_right) {
        datR <- BayesGLM_results$cortex_right$result_multiple[[meas]]*1
        mwallR <- !is.na(datR)
        datR <- datR[mwallR]
      }
      if (do_sub) {
        datSub <- BayesGLM_results$subcortical$result_multiple[[meas]]*1
        colnames(datSub) <- NULL
      }
      xii_meas <- as.xifti(
        cortexL = datL,
        cortexL_mwall = mwallL,
        cortexR = datR,
        cortexR_mwall = mwallR,
        subcortVol = datSub,
        subcortLabs = submeta_ss$labels,
        subcortMask = submeta_ss$mask
      )
      xii_meas$meta$subcort$trans_mat <- submeta_ss$trans_mat
      xii_meas$meta$subcort$trans_units <- submeta_ss$trans_units
      xii_meas$meta$cifti$names <- field_names
      result_multiple_xii[[meas]] <- xii_meas
    }

    #stuff we don't have when fitting multiple models
    HRFs <- FIR <- design_FIR <- stimulus <- NULL
    #field_names <- NULL
  }

  result <- list(
    field_estimates_xii = list(
      Bayes = field_cifti,
      classical = field_cifti_classical
    ),
    result_multiple_xii = result_multiple_xii,
    design = design, # after centering/scaling, before nuisance regression / prewhitening
    design_multiple = design_multiple,
    nuisance = nuisance,
    stimulus = stimulus,
    HRFs = HRFs, #convolved HRFs
    HRF_basis = HRF, #HRF basis functions
    FIR = FIR,
    design_FIR = design_FIR,
    field_names = field_names,
    session_names = session_names,
    n_sess = nS_orig,
    BayesGLM_results = BayesGLM_results
  )

  if (verbose>0) cat('Done!\n')

  class(result) <- "BayesGLM_cifti"
  result
}
