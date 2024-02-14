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
#' @param design See \code{make_design}.
#' @param dHRF \code{dHRF} controls the addition of HRF derivatives to the design matrix if
#'  \code{onsets} is provided. Set to \code{0} to not model the derivatives and
#'  only include the main HRF regressor; set to \code{1} (default) to model the
#'  temporal derivative too; or, set to \code{2} to model both the temporal
#'  derivative and the dispersion derivative too.
#'  If \code{dHRF==0}, there is one design column per task. If \code{dHRF==1},
#'  there are two design columns per task. And if \code{dHRF==2}, there are
#'  three design columns per task. If there are several tasks and \code{dHRF>0},
#'  spatial modeling with INLA may require large computation times. A possible
#'  adjustment is to model the columns for HRF derivatives as nuisance signals
#'  rather than fields. This can be controlled by the \code{dHRF_as} argument
#'  to \code{BayesGLM(_cifti)}.
#' @param dHRF_as \code{"auto"} (default) ...
#' @inheritParams session_names_Param
#' @param nuisance (Optional) A \eqn{T \times J} matrix of nuisance signals.
#'  These are regressed from the fMRI data and the design matrix prior to the
#'  GLM computation. For multi-session modeling, this argument should be a list
#'  of such matrices.
#' @param TR The temporal resolution of the data, in seconds.
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
#' @importFrom abind abind
#'
#' @export
BayesGLM_cifti <- function(
  cifti_fname, design,
  surfL_fname=NULL,
  surfR_fname=NULL,
  brainstructures=c('left','right'),
  dHRF=c(1, 0, 2),
  dHRF_as=c("auto", "nuisance", "field"),
  session_names = NULL,
  nuisance=NULL,
  TR=NULL, hpf=NULL,
  DCT=if(is.null(hpf)) {4} else {NULL},
  resamp_res=10000,
  nbhd_order=1, buffer=c(1,1,3,4,4),
  # Below arguments shared with `BayesGLM`.
  #combine_sessions = TRUE,
  scale_BOLD = c("auto", "mean", "sd", "none"),
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

  dHRF <- as.numeric(match.arg(as.character(dHRF), c("1", "0", "2")))
  dHRF_as <- match.arg(dHRF_as, c("auto", "nuisance", "field"))

  # Nuisance arguments.
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
    if (!is.null(nuisance)) nuisance <- list(single_session = nuisance)

  } else {

    if (verbose>0) { cat(paste0(
      "Preparing to analyze ",nS,
      " task fMRI sessions with a common set of tasks.\n"
    )) }

    #name sessions
    if (is.null(session_names)) session_names <- paste0('session', 1:nS)

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

  # Data setup. ----------------------------------------------------------------
  if (verbose>0) cat('Setting up data:\n')

  cbind2 <- function(mat_or_NULL, to_add) {
    if (!is.null(mat_or_NULL)) {
      cbind(mat_or_NULL, to_add)
    } else {
      to_add
    }
  }

  ## xifti things. -------------------------------------------------------------
  ### For each session, separate the CIFTI data into left/right/sub and read in files
  BOLD_list <- list(left=NULL, right=NULL, subcortical=NULL)
  mwallL <- mwallR <- NULL
  nT <- vector("numeric", nS)

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
    nT[ss] <- ncol(xii_ss)

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

  ## Design and HRF derivatives. -----------------------------------------------
  stopifnot(is.list(design))
  stopifnot(is.list(design$design))
  if (length(design$design) != nS) { stop("The length of `design` should match the number of sessions, ", nS, ".") }
  design_is_multiple <- length(dim(design$design[[1]])) == 3
  nD <- if (design_is_multiple) { dim(design$design[[1]])[3] } else { 1 }
  design_is_from_onsets <- !is.null(design[[1]]$onsets_misc)
  # [TO DO]: more rigorous checks?

  # Determine what to do with the HRF derivatives, if `dHRF_as=="auto`
  # Determine the number of design matrix columns.
  if (design_is_from_onsets && dHRF > 0) {
    nJ <- dim(design$design[[1]])[2]
    nK_if_add_derivs <- (nJ * (seq(0,2) + 1))[dHRF + 1]
    if (dHRF_as=="auto") {
      dHRF_as <- ifelse(nK_if_add_derivs > 5, "nuisance", "field")
      message("Using the HRF derivatives as ", switch(dHRF_as,
        nuisance="nuisance signals to remove.",
        field="fields in the model."
      ))
    }
    nK <- switch(dHRF_as, nuisance=nJ, field=nK_if_add_derivs)
  } else {
    nK <- dim(design$design[[1]])[2]
  }

  # Warn the user if the number of design matrix columns exceeds five.
  if (Bayes && nK > 5) {
    message("The number of regressors to be modeled spatially exceeds five. INLA computation may be slow. Consider reducing the number of design matrix columns, e.g. by modeling HRF derivatives as nuisance")
    Sys.sleep(10)
  }

  # Move the HRF derivatives to the design or nuisance matrices, if applicable.
  if (dHRF > 0) {
    for (ss in seq(nS)) {
      if (dHRF_as=="design") {
        to_add <- switch(dHRF,
          design[[ss]]$onsets_misc$HRFs_d,
          cbind(design[[ss]]$onsets_misc$HRFs_d, design[[ss]]$onsets_misc$HRFs_dd)
        )
        if (design_is_multiple) {
          to_add <- do.call(cbind, rep(list(to_add), nD))
          design[[ss]]$design <- abind::abind(design, to_add, along=1)
        } else {
          design[[ss]]$design <- cbind2(design[[ss]]$design, to_add)
        }
      } else if (dHRF_as=="nuisance") {
        to_add <- switch(dHRF,
          design[[ss]]$onsets_misc$HRFs_d,
          cbind(design[[ss]]$onsets_misc$HRFs_d, design[[ss]]$onsets_misc$HRFs_dd)
        )
        nuisance[[ss]] <- cbind2(nuisance[[ss]], to_add)
      }
    }
  }

  # Get `field_names` and final design matrix. Separate out misc. onsets info.
  field_names <- colnames(design$design[[1]])
  if (design_is_from_onsets) {
    stimulus <- design$onsets_misc$stimulus
    HRFs <- design$onsets_misc$HRFs #convolved HRFs
    HRF_basis <- design$onsets_misc$HRF #HRF basis functions
    FIR <- design$onsets_misc$FIR
    design_FIR <- design$onsets_misc$design_FIR
  } else {
    stimulus <- HRFs <- HRF_basis <- FIR <- design_FIR <- NULL
  }
  design_onsets_misc <- design$onsets_misc[c("HRFs", "HRFs_d", "HRFs_dd")]
  design_onsets_misc <- lapply(seq(nS), function(ss){lapply(design_onsets_misc, '[[', ss)})
  design <- design$design

  ## Nuisance matrices. --------------------------------------------------------
  # Add DCT bases to nuisance matrix
  DCTs <- vector("numeric", nS)
  for (ss in 1:nS) {
    # DCT highpass filter
    if (!is.null(hpf) || !is.null(DCT)) {
      # Get the num. of bases for this session.
      if (!is.null(hpf)) {
        DCTs[ss] <- round(dct_convert(nT[ss], TR, f=hpf))
        if(verbose > 0 & ss==1) cat(paste('\tIncluding',DCTs[ss],'DCT basis functions in the model based for hpf =',hpf,'Hz\n'))
      } else {
        DCTs[ss] <- DCT
      }
      # Generate the bases and add them.
      DCTb_ss <- dct_bases(nT[ss], DCTs[ss])
      if (DCTs[ss] > 0) { nuisance[[ss]] <- cbind2(nuisance[[ss]], DCTb_ss) }
    }
  }

  # Do GLM. --------------------------------------------------------------------
  BayesGLM_results <- list(left = NULL, right = NULL, subcortical = NULL)

  ## Loop through brainstructures. ---------------------------------------------
  for (bb in brainstructures) {
    bname <- list(left="left cortex", right="right cortex", subcortical="subcortical")[[bb]]
    if (verbose>0) cat(paste0(toupper(bname)," analysis:\n"))

    # Set up session list.
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

    # Set up spatial information.
    vertices <- faces <- labels <- NULL
    if(do_cort){
      vertices <- spatial_list[[bb]]$vertices
      faces <- spatial_list[[bb]]$faces
    } else {
      labels <- spatial_list[[bb]]
    }

    # `BayesGLM` call.
    BayesGLM_out <- BayesGLM(
      data = session_data,
      #design_multiple = design_multiple,
      vertices = vertices,
      faces = faces,
      labels = labels,
      nbhd_order = nbhd_order,
      buffer = buffer,
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

    # # Update session info if averaged over sessions.
    # if (bb == brainstructures[1] && combine_sessions) {
    #   session_names <- BayesGLM_out$session_names
    #   nS <- 1
    # }

    # Clear up memory for next iteration.
    rm(BayesGLM_out); gc()
  }
  names(BayesGLM_results)[names(BayesGLM_results)=="left"] <- "cortex_left"
  names(BayesGLM_results)[names(BayesGLM_results)=="right"] <- "cortex_right"

  ## Construct beta estimates as `xifti` objects. ------------------------------

  if (verbose>0) cat("Formatting results.\n")

  results_xii <- list(classical=vector('list', nS), Bayesian=vector('list', nS))
  bestmodel_xii <- sigma2_xii <- vector('list', nS)
  names(results_xii$classical) <- names(results_xii$Bayesian) <- names(bestmodel_xii) <- names(sigma2_xii) <- session_names

  if (!design_is_multiple) {
    datL <- datR <- datSub <- NULL
    for (method in c("classical", "Bayesian")[seq(1+do_Bayesian)]) {
      for (ss in seq(nS)) {
        if (do_left) {
          datL <- switch(method,
            classical = BayesGLM_results$cortex_left$result_classical[[ss]]$estimates,
            Bayesian = BayesGLM_results$cortex_left$field_estimates[[ss]]
          )
          # Update mwall b/c `mask2` in `BayesGLM` can change the medial wall.
          mwallL <- !is.na(datL[,1])
          datL <- datL[mwallL,]
          colnames(datL) <- NULL
        }
        if (do_right) {
          datR <- switch(method,
            classical = BayesGLM_results$cortex_right$result_classical[[ss]]$estimates,
            Bayesian = BayesGLM_results$cortex_right$field_estimates[[ss]]
          )
          mwallR <- !is.na(datR[,1])
          datR <- datR[mwallR,]
          colnames(datR) <- NULL
        }
        if (do_sub) {
          datSub <- switch(method,
            classical = BayesGLM_results$subcort$result_classical[[ss]]$estimates,
            Bayesian = BayesGLM_results$subcortical$field_estimates[[ss]]
          )
          colnames(datSub) <- NULL
        }
        results_xii[[method]][[ss]] <- as.xifti(
          cortexL = datL,
          cortexL_mwall = if (do_left) { mwallL } else { NULL },
          cortexR = datR,
          cortexR_mwall = if (do_right) { mwallR } else { NULL },
          subcortVol = datSub,
          subcortLabs = submeta_ss$labels,
          subcortMask = submeta_ss$mask
        )
        results_xii[[method]][[ss]]$meta$subcort$trans_mat <- submeta_ss$trans_mat
        results_xii[[method]][[ss]]$meta$subcort$trans_units <- submeta_ss$trans_units
        results_xii[[method]][[ss]]$meta$cifti$names <- field_names
      }
    }

    names(design) <- names(nuisance) <- session_names
    if (!is.null(stimulus)) { names(stimulus) <- session_names }
    #if(length(FIR) > 0) names(FIR) <- names(design_FIR) <- session_names #this check doesn't work for multi-session modeling

  } else {
    datL <- datR <- datSub <- NULL #index of best-fitting model
    betaL <- betaR <- betaSub <- NULL #beta estimates for best-fitting model
    sigma2L <- sigma2R <- sigma2Sub <- NULL #residual var of models

    for (ss in seq(nS)) {
      # INDEX OF BEST MODEL
      if (do_left) {
        datL <- BayesGLM_results$cortex_left$result_multiple[[ss]]$bestmodel #index of best model
        betaL <- BayesGLM_results$cortex_left$result_multiple[[ss]]$beta_estimates #V x K x P (P = number of models tested)
        sigma2L <- BayesGLM_results$cortex_left$result_multiple[[ss]]$sigma2 #V x P (P = number of models tested)
        #only save beta estimates for the best fitting model
        betaL <- apply(betaL, 1, as.matrix, simplify=FALSE) #form into a list of length V, each a K x P matrix
        betaL <- t(mapply(function(matrix, index) matrix[, index, drop = FALSE], betaL, datL, SIMPLIFY = TRUE)) #beta estimates (VxK) for the best model
        # Update mwall b/c `mask2` in `BayesGLM` can change the medial wall.
        mwallL <- !is.na(datL)
        datL <- datL[mwallL]
        betaL <- betaL[mwallL,]
        sigma2L <- sigma2L[mwallL,]
      }
      if (do_right) {
        datR <- BayesGLM_results$cortex_right$result_multiple[[ss]]$bestmodel
        betaR <- BayesGLM_results$cortex_right$result_multiple[[ss]]$beta_estimates #V x K x P (P = number of models tested)
        sigma2R <- BayesGLM_results$cortex_right$result_multiple[[ss]]$sigma2 #V x P (P = number of models tested)
        #only save beta estimates for the best fitting model
        betaR <- apply(betaR, 1, as.matrix, simplify=FALSE) #form into a list of length V, each a K x P matrix
        betaR <- t(mapply(function(matrix, index) matrix[, index, drop = FALSE], betaR, datR, SIMPLIFY = TRUE)) #beta estimates for the best model
        mwallR <- !is.na(datR)
        datR <- datR[mwallR]
        betaR <- betaR[mwallR,]
        sigma2R <- sigma2R[mwallR,]
      }
      if (do_sub) {
        #[TO DO]: do this for subcortex, as for L and R above
        #datSub <- BayesGLM_results$subcortical$result_multiple[[ss]]$bestmodel
        #colnames(datSub) <- NULL
      }

      bestmodel_xii[[ss]] <- as.xifti(
        cortexL = datL,
        cortexL_mwall = if (do_left) { mwallL } else { NULL },
        cortexR = datR,
        cortexR_mwall = if (do_right) { mwallR } else { NULL },
        subcortVol = datSub,
        subcortLabs = submeta_ss$labels,
        subcortMask = submeta_ss$mask
      )
      bestmodel_xii[[ss]]$meta$subcort$trans_mat <- submeta_ss$trans_mat
      bestmodel_xii[[ss]]$meta$subcort$trans_units <- submeta_ss$trans_units
      bestmodel_xii[[ss]]$meta$cifti$names <- field_names

      results_xii$classical[[ss]] <- as.xifti(
        cortexL = betaL,
        cortexL_mwall = if (do_left) { mwallL } else { NULL },
        cortexR = betaR,
        cortexR_mwall = if (do_right) { mwallR } else { NULL },
        subcortVol = betaSub,
        subcortLabs = submeta_ss$labels,
        subcortMask = submeta_ss$mask
      )
      results_xii$classical[[ss]]$meta$subcort$trans_mat <- submeta_ss$trans_mat
      results_xii$classical[[ss]]$meta$subcort$trans_units <- submeta_ss$trans_units
      results_xii$classical[[ss]]$meta$cifti$names <- field_names

      sigma2_xii[[ss]] <- as.xifti(
        cortexL = sigma2L,
        cortexL_mwall = if (do_left) { mwallL } else { NULL },
        cortexR = sigma2R,
        cortexR_mwall = if (do_right) { mwallR } else { NULL },
        subcortVol = sigma2Sub,
        subcortLabs = submeta_ss$labels,
        subcortMask = submeta_ss$mask
      )
      sigma2_xii[[ss]]$meta$subcort$trans_mat <- submeta_ss$trans_mat
      sigma2_xii[[ss]]$meta$subcort$trans_units <- submeta_ss$trans_units
      sigma2_xii[[ss]]$meta$cifti$names <- field_names
    }

    #stuff we don't have when fitting multiple models
    HRFs <- FIR <- design_FIR <- stimulus <- NULL
    field_names <- field_names
  }

  result <- list(
    estimates_xii = list(
      Bayes = results_xii$Bayesian,
      classical = results_xii$classical
    ),
    bestmodel_xii = bestmodel_xii, #only has values for multiple design comparison case
    sigma2_xii = sigma2_xii, #only has values for multiple design comparison case
    design = design, # after centering/scaling, before nuisance regression / prewhitening
    #design_multiple = design_multiple,
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
