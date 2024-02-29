#' BayesGLM for CIFTI
#'
#' Performs spatial Bayesian GLM for task fMRI activation with CIFTI-format
#'  data. The cortex is modeled as a surface mesh, and subcortical structures
#'  are modeled as distinct volumetric regions. Includes the pre-processing
#'  steps of nuisance regression, prewhitening, scaling, and variance
#'  normalization. Supports both single- and multi-session analysis. Can also
#'  compute just the classical (spatially-independent) GLM. Can also compare
#'  different choices of designs for single-session classical GLM.
#'
#' To use \code{BayesGLM_cifti}, the design matrix must first be constructed
#'  with \code{\link{make_design}}.
#'
#' @inheritSection Connectome_Workbench_Description Connectome Workbench Requirement
#' @inheritSection INLA_Description INLA Requirement
#' @inheritSection INLA_Latent_Fields_Limit_Description INLA Latent Fields Limit
#'
#' @inheritParams BOLD_param_BayesGLM_cifti
#' @param design A \code{"BayesfMRI_design"} object from \code{\link{make_design}}.
#' @inheritParams brainstructures_param_BayesGLM_cifti
#' @inheritParams surfaces_Param_BayesGLM_cifti
#' @inheritParams resamp_res_Param_BayesGLM_cifti
#' @inheritParams nbhd_order_Param
#' @inheritParams buffer_Param
#' @inheritParams nuisance_param_BayesGLM_cifti
#' @inheritParams detrending_param_BayesGLM_cifti
#' @inheritParams scale_BOLD_Param
#' @inheritParams Bayes_Param
# @inheritParams EM_Param
#' @inheritParams ar_order_Param
#' @inheritParams ar_smooth_Param
#' @inheritParams aic_Param
#' @inheritParams n_threads_Param
#' @inheritParams return_INLA_Param
#' @inheritParams verbose_Param
#' @inheritParams mean_var_Tol_Param
# @inheritParams emTol_Param
#'
#' @return An object of class \code{"BayesGLM_cifti"}: a list with elements
#'  \describe{
#'    \item{betas_Bayesian}{The field coefficients for the Bayesian model.}
#'    \item{betas_classical}{The field coefficients for the classical model.}
#'    \item{GLMs_Bayesian}{The entire list of GLM results, except for parameters estimated for the classical model.}
#'    \item{GLMs_classical}{Parameters estimated for the classical model from the GLM.}
#'    \item{brainstructures}{\code{data.frame} summarizing the spatial features of each brain structure modeled.}
#'    \item{sessions}{\code{data.frame} with the \code{name} and \code{nTime} of each BOLD session.}
#'    \item{fields}{\code{data.frame} with the \code{name}, related \code{task}, and \code{HRF_order} of each field.}
#'  }
#'
# @importFrom ciftiTools read_cifti resample_gifti as.xifti remove_xifti
#' @import ciftiTools
#' @importFrom fMRItools unmask_mat match_input
#' @importFrom matrixStats rowVars rowSums2 colVars
#' @importFrom parallel detectCores
#' @importFrom Matrix bdiag
#'
#' @export
#'
BayesGLM_cifti <- function(
  BOLD,
  design,
  brainstructures=c('left','right'),
  surfL=NULL,
  surfR=NULL,
  resamp_res=10000,
  nbhd_order=1,
  buffer=c(1,1,3,4,4),
  nuisance=NULL,
  TR=NULL, hpf=NULL,
  DCT=if(is.null(hpf)) {4} else {NULL},
  # Below arguments shared with `BayesGLM`.
  scale_BOLD = c("auto", "mean", "sd", "none"),
  Bayes = TRUE,
  #EM = FALSE,
  ar_order = 6,
  ar_smooth = 5,
  aic = FALSE,
  n_threads = 4,
  return_INLA = c("trimmed", "full", "minimal"),
  verbose = 1,
  meanTol = 1e-6,
  varTol = 1e-6#,emTol = 1e-3
){

  EM <- FALSE
  emTol <- 1e-3

  # Preliminaries. -------------------------------------------------------------
  ### Simple parameter checks. -------------------------------------------------

  # `design`: check it's from `make_design` and grab the data dimensions.
  stopifnot(inherits(design, "BfMRI_design"))
  des_type <- design$design_type
  if (des_type == "compare") { stop("Use `CompareGLM_cifti` to compare models.") }
  des_dims <- setNames(design$dims$count, rownames(design$dims))
  nS <- as.numeric(des_dims["sessions"])
  session_names <- names(design$design)
  nT <- vapply(design$design, function(q){dim(q)[1]}, 0)
  nK <- as.numeric(des_dims["fields"])
  field_names <- design$field_names
  if (Bayes && nK > 5) {
    message(
      "The number of regressors to be modeled spatially exceeds five. ",
      "INLA computation may be slow. Consider reducing the number of design ",
      "matrix columns, e.g. by modeling HRF derivatives as nuisance. See ",
      "`remove_from_design` for a helper function for this task."
    )
    Sys.sleep(10)
  }
  nD <- as.numeric(des_dims["design_matrices"])

  stopifnot(is.numeric(nbhd_order))
  stopifnot(fMRItools::is_1(nbhd_order, "numeric"))
  stopifnot(nbhd_order>0 && nbhd_order==round(nbhd_order))
  stopifnot(is.numeric(buffer))

  do <- vector("list")

  # In a separate function because these checks are shared with `BayesGLM`.
  x <- BayesGLM_argChecks(
    scale_BOLD = scale_BOLD,
    Bayes = Bayes,
    EM = EM,
    ar_order = ar_order,
    ar_smooth = ar_smooth,
    aic = aic,
    n_threads = n_threads,
    return_INLA = return_INLA,
    verbose = verbose,
    meanTol = meanTol,
    varTol = varTol,
    emTol = emTol
  )
  scale_BOLD <- x$scale_BOLD
  do$Bayesian <- x$Bayes; rm(Bayes) # rename
  do$EM <- x$do_EM; rm(EM) # rename
  do$pw <- x$do_pw # unused
  return_INLA <- x$return_INLA
  rm(x)

  ### Brain structures. --------------------------------------------------------
  if ("both" %in% brainstructures) { brainstructures <- c("left", "right") }
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }
  brainstructures <- fMRItools::match_input(
    brainstructures, c("left","right","subcortical"),
    user_value_label="brainstructures"
  )
  do$left <- ('left' %in% brainstructures)
  do$right <- ('right' %in% brainstructures)
  do$sub <- ('subcortical' %in% brainstructures)
  do$cortex <- do$left || do$right
  if (!do$cortex) { resamp_res <- NULL }

  ### Initialize `spatial` to store all spatial information.
  spatial <- list(
    cortexL = list(surf=NULL, mask=NULL),
    cortexR = list(surf=NULL, mask=NULL),
    subcort = list(
      label=NULL,
      trans_mat=NULL, trans_units=NULL,
      nbhd_order=nbhd_order, buffer=buffer
    )
  )
  if (!do$left) { spatial$cortexL <- NULL }
  if (!do$right) { spatial$cortexR <- NULL }
  if (!do$sub) { spatial$subcort <- NULL }

  ### Initialize `nV0` to store the total number of locations for each brain structure.
  nV0 <- as.character(lapply(spatial, function(q){NA}))

  ### Checks for `BOLD` that don't require reading the data in. ----------------
  # Make `BOLD` a sessions-length character vector, or a sessions-length list of
  #  \code{"xifti"} objects.
  is_xifti <- FALSE
  if (is.character(BOLD)) {
    BOLD <- as.list(BOLD)
  } else if (is.xifti(BOLD, messages=FALSE)) {
    is_xifti <- TRUE
    BOLD <- list(BOLD)
  } else if (is.list(BOLD)) {
    if (all(vapply(BOLD, is.character, FALSE)) && all(vapply(BOLD, length, 0)==1)) {
      BOLD <- setNames(as.character(BOLD), names(BOLD))
    } else {
      is_xifti_vec <- vapply(BOLD, is.xifti, messages=FALSE, FALSE)
      if (!all(is_xifti_vec)) {
        stop('`BOLD` should be a character vector or list of `"xifti"` objects.')
      }
      rm(is_xifti_vec)
      is_xifti <- TRUE
    }
  } else {
    stop('`BOLD` should be a character vector or list of `"xifti"` objects.')
  }

  if (length(BOLD) != nS) {
    stop(
      "The length of `BOLD`, ", length(BOLD),
      " does not match the number of sessions indicated by `design`, ", nS, "."
    )
  }
  if (!is.null(names(BOLD)) && !all(names(BOLD) == session_names)) {
    warning("Using session names from `design`; ignoring BOLD names.")
  }
  names(BOLD) <- session_names

  if (verbose > 0) {
    if (nS==1) {
      cat("Preparing to analyze a single task fMRI session.\n")
    } else {
      cat("Preparing to analyze", nS, "task fMRI sessions with a common set of tasks.\n")
    }
  }

  ### Surfaces for cortex analysis. --------------------------------------------
  # Use surfaces in this order of priority: `surfL` parameter; surface included
  #   in `BOLD` `xifti`; HCP group average inflated surface from `ciftiTools`.
  # Still do this even if `!Bayes` to get surfaces for visualization.
  # Note that surfaces may be resampled (again)
  if (do$left) {
    # Get surface.
    if (is.null(surfL)) {
      if (is_xifti && !is.null(BOLD[[1]]$surf$cortex_left)) {
        all_left_surfs <- c(
          lapply(BOLD, function(q){q$surf$cortex_left}), list(NULL)
        )
        if (length(unique(all_left_surfs)) > 2) {
          warning("Using left surface from the first `BOLD` `xifti` for all modeling. Ignoring the other left surfaces.\n")
        } else {
          cat("Using left surface from the `BOLD` `xifti` data.\n")
        }
        surfL <- BOLD[[1]]$surf$cortex_left
      } else {
        surfL <- ciftiTools.files()$surf["left"]
      }
    }
    # Read and resample, if necessary.
    if (suppressMessages(is.surf(surfL))) {
      if (!is.null(resamp_res)) { surfL <- resample_surf(surfL, resamp_res=resamp_res) }
    } else {
      surfL <- read_surf(surfL, resamp_res=resamp_res)
    }
    spatial$cortexL$surf <- surfL; rm(surfL)
  }

  if (do$right) {
    # Get surface.
    if (is.null(surfR)) {
      if (is_xifti && !is.null(BOLD[[1]]$surf$cortex_right)) {
        all_right_surfs <- c(
          lapply(BOLD, function(q){q$surf$cortex_right}), list(NULL)
        )
        if (length(unique(all_right_surfs)) > 2) {
          warning("Using right surface from the first `BOLD` `xifti` for all modeling. Ignoring the other right surfaces.\n")
        } else {
          cat("Using right surface from the `BOLD` `xifti` data.\n")
        }
        surfR <- BOLD[[1]]$surf$cortex_right
      } else {
        surfR <- ciftiTools.files()$surf["right"]
      }
    }
    # Read and resample, if necessary.
    if (suppressMessages(is.surf(surfR))) {
      if (!is.null(resamp_res)) { surfR <- resample_surf(surfR, resamp_res=resamp_res) }
    } else {
      surfR <- read_surf(surfR, resamp_res=resamp_res)
    }
    spatial$cortexR$surf <- surfR; rm(surfR)
  }

  # `BOLD`: read in data and metadata, for each brain structure. ---------------

  BOLD_input_msg <- function(ss, nS, do=c("read", "resample")){
    do <- switch(do, read="Reading", resample="Resampling")
    out <- if (ss==1) { paste0("\t", do,  " BOLD data") } else { "" }
    if (nS==1) {
      out <- paste0(out, ".\n")
    } else {
      if (ss==1) { out <- paste0(out, " for session ") }
      out <- paste0(out, ss, ifelse(ss!=nS, ", ", ".\n"))
    }
    out
  }

  for (ss in seq(nS)) {
    ### Read and/or resample the CIFTI data. -----------------------------------
    if (is_xifti) {
      if (do$cortex) {
        if (!is.null(resamp_res)) {
          if (any(ciftiTools::infer_resolution(BOLD[[ss]])!=resamp_res)) {
            if (verbose>0) { cat(BOLD_input_msg(ss, nS, "resample")) }
            BOLD[[ss]] <- resample_xifti(BOLD[[ss]], resamp_res=resamp_res)
          }
        }
      }
    } else {
      if (verbose>0) { cat(BOLD_input_msg(ss, nS, "read")) }
      BOLD[[ss]] <- read_cifti(
        BOLD[[ss]], brainstructures=brainstructures,
        resamp_res=resamp_res
      )
    }

    if (ncol(BOLD[[ss]]) != nT[ss]) { stop(
      "The design for session ", session_names[ss], " indicates ", nT[ss], " ",
      "volumes, but the `xifti` data for this session has ", ncol(BOLD[[ss]]), " ",
      "volumes. Repeat `make_design` with a revised `nTime`, or correct `BOLD`."
    )}

    if (do$left && is.null(BOLD[[ss]]$data$cortex_left)) { stop("Left cortex data is missing from this BOLD session.") }
    if (do$right && is.null(BOLD[[ss]]$data$cortex_right)) { stop("Right cortex data is missing from this BOLD session.") }
    if (do$sub && is.null(BOLD[[ss]]$data$subcort)) { stop("Subcortex data is missing from this BOLD session.") }

    ### Total spatial dims, `nV0`. Check it's consistent across sessions. ------
    # `xii_res`: total spatial dims according to the BOLD data.
    #   Names: 'left', 'right', 'subcort'.
    xii_res <- vector("numeric", 0)
    if (do$cortex) {
      xii_res <- c(xii_res, ciftiTools::infer_resolution(BOLD[[ss]]))
    }
    if (do$sub) {
      xii_res <- c(xii_res, c(subcort=sum(BOLD[[ss]]$meta$subcort$mask)))
    }

    # Set `nV0` based on the `xii_res` of the first session...
    if (ss == 1) {
      if (do$left) {
        # Check that the left cortex data is nonempty.
        if (xii_res["left"] %in% c(0, NA, NaN)) { stop("This BOLD session does not seem to have left cortex data.") }
        surfL_res <- get_nV(spatial$cortexL)$T
        # Check that the left cortex resolution matches the surface resolution.
        if (xii_res["left"] != get_nV(spatial$cortexL)$T) {
          if (is.null(resamp_res)) {
            # If mismatch and `resamp_res` was `NULL`, try resampling the surface to match the data.
            spatial$cortexL$surf <- resample_surf(spatial$cortexL$surf, xii_res["left"])
            if(xii_res["left"] != get_nV(spatial$cortexL)$T) {
              stop("The left surface could not be resampled to the resolution of the left cortex BOLD data,", xii_res["left"])
            }
          } else {
            # If mismatch and `resamp_res` was not `NULL`, we have a problem.
            stop("The left surface could not be resampled to the resolution of the left cortex BOLD data,", xii_res["left"])
          }
        }
        nV0["cortexL"] <- xii_res["left"]
      }
      if (do$right) {
        # Check that the right cortex data is nonempty.
        if (xii_res["right"] %in% c(0, NA, NaN)) { stop("This BOLD session does not seem to have right cortex data.") }
        surfL_res <- get_nV(spatial$cortexR)$T
        # Check that the right cortex resolution matches the surface resolution.
        if (xii_res["right"] != get_nV(spatial$cortexR)$T) {
          if (is.null(resamp_res)) {
            # If mismatch and `resamp_res` was `NULL`, try resampling the surface to match the data.
            spatial$cortexR$surf <- resample_surf(spatial$cortexR$surf, xii_res["right"])
            if(xii_res["right"] != get_nV(spatial$cortexR)$T) {
              stop("The right surface could not be resampled to the resolution of the right cortex BOLD data,", xii_res["right"])
            }
          } else {
            # If mismatch and `resamp_res` was not `NULL`, we have a problem.
            stop("The right surface could not be resampled to the resolution of the right cortex BOLD data,", xii_res["left"])
          }
        }
        nV0["cortexR"] <- xii_res["right"]
      }
      if (do$sub) {
        nV0["subcort"] <- sum(BOLD[[ss]]$meta$subcort$mask)
      }

      # Per-location design: check `sum(nV0)` matches with `design`.
      if (des_type == "per_loc") {
        if (sum(nV0) != dim(design$design[[1]])[3]) { stop(
          "`design` indicates ", dim(design$design[[1]])[3], " ",
          "total locations, but the `xifti` data for this session has ",
          sum(nV0), " total locations. Repeat `make_design` with a corrected ",
          "design, or fix `BOLD`."
        )}
      }

    # Check `nV0` matches `xii_res` of other sessions.
    } else {
      if (do$left && (xii_res["left"] != nV0["cortexL"])) {
        stop(
          "This BOLD session appears to have left cortex data with ",
          xii_res["left"], "total vertices, while the first session ",
          "appears to have", nV0["cortexL"], "total vertices."
        )
      }
      if (do$right && (xii_res["right"] != nV0["cortexR"])) {
        stop(
          "This BOLD session appears to have right cortex data with ",
          xii_res["right"], "total vertices, while the first session ",
          "appears to have", nV0["cortexR"], "total vertices."
        )
      }
      if (do$sub && (xii_res["subcort"] != nV0["subcort"])) {
        stop(
          "This BOLD session appears to have subcortical data with",
          xii_res["subcort"], "total voxels, while the first session ",
          "appears to have", nV0["subcort"], " total voxels."
        )
      }
    }
    rm(xii_res)

    ### Collect `spatial` metadata. --------------------------------------------
    # Separate BOLD data and metadata.
    xii_meta <- BOLD[[ss]]$meta
    BOLD[[ss]] <- BOLD[[ss]]$data

    # cortex: `mask`. subcortex: `label`, `trans_mat`, `trans_units`.
    if (ss == 1) {
      if (do$left) {
        spatial$cortexL$mask <- xii_meta$cortex$medial_wall_mask$left
        if (is.null(spatial$cortexL$mask)) { spatial$cortexL$mask <- rep(TRUE, nV0["cortexL"]) }
      }
      if (do$right) {
        spatial$cortexR$mask <- xii_meta$cortex$medial_wall_mask$right
        if (is.null(spatial$cortexR$mask)) { spatial$cortexR$mask <- rep(TRUE, nV0["cortexR"]) }
      }
      if (do$sub) {
        submeta <- xii_meta$subcort
        spatial$subcort$label <- submeta$mask*0;
        spatial$subcort$label[submeta$mask==TRUE] <- submeta$labels
        spatial$subcort$trans_mat <- submeta$trans_mat
        spatial$subcort$trans_units <- submeta$trans_units
      } else {
        submeta <- NULL
      }
    # For multi-session: ensure medial walls and subcortex masks match.
    # [TO DO]: just use their intersection.
    } else {
      if (do$left) {
        stopifnot(all(spatial$cortexL$mask == xii_meta$cortex$medial_wall_mask$left))
      }
      if (do$right) {
        stopifnot(all(spatial$cortexR$mask == xii_meta$cortex$medial_wall_mask$right))
      }
      if (do$sub) {
        stopifnot(length(dim(spatial$subcort$label)) == length(dim(xii_meta$subcort$mask)))
        stopifnot(all(dim(spatial$subcort$label)) == dim(xii_meta$subcort$mask))
        stopifnot(all((spatial$subcort$label!=0) == xii_meta$subcort$mask))
      }
    }
    rm(xii_meta)
  }

  ### Collate `BOLD` by brainstructure. ----------------------------------------
  BOLD <- lapply(
    c("cortex_left", "cortex_right", "subcort"),
    function(bs){ setNames(lapply(BOLD, '[[', bs), session_names) }
  )
  names(BOLD) <- c("cortexL", "cortexR", "subcort")
  if (!do$left) { BOLD$cortexL <- NULL }
  if (!do$right) { BOLD$cortexR <- NULL }
  if (!do$sub) { BOLD$subcort <- NULL }

  # Transpose all `BOLD` to `TxV`.
  BOLD <- lapply(BOLD, function(q){lapply(q, t)})

  ### Check `nuisance`. --------------------------------------------------------
  if (!is.null(nuisance)) {
    # Make `nuisance` an `nS`-length list.
    if (!is.list(nuisance)) {
      if (nS==1) { nuisance = list(nuisance) } else { stop(
        "If provided, `nuisance` should be a sessions-length list."
      )}
    } else {
      if (length(nuisance) != nS) {
        "If provided, `nuisance` should be a sessions-length list."
      }
      if (!is.null(names(nuisance)) && !all(names(nuisance) == session_names)) {
        #warning("Ignoring `names(nuisance)`; use `session_names` in `make_design`.")
      }
    }
    names(nuisance) <- session_names
    # Ensure each element is a numeric matrix with `nT` rows.
    for (ss in seq(length(nS))) {
      stopifnot(is.matrix.or.df(nuisance[[ss]]))
      stopifnot(is.numeric(nuisance[[ss]]))
      stopifnot(nrow(nuisance[[ss]]) == nT[ss])
    }
  } else {
    nuisance <- vector("list", nS)
  }

  ### Make DCT bases in `design` for the high-pass filter. ---------------------
  for (ss in seq(nS)) {
    nuisance[[ss]] <- cbind2(nuisance[[ss]],
      BayesGLM_cifti_make_DCT(DCT, hpf, nT[ss], verbose)
    )
  }

  # Do GLM. --------------------------------------------------------------------
  BayesGLM_results <- list(cortexL = NULL, cortexR = NULL, subcort = NULL)

  ## Loop through brainstructures. ---------------------------------------------
  bs_names <- data.frame(
    d = c("cortexL", "cortexR", "subcort"), # Names used in this code.
    v = c("Left cortex", "Right cortex", "Subcortex") # Verbose names to show user.
  )

  for (bb in seq(nrow(bs_names))) {
    if (!(bs_names$d[bb] %in% names(BOLD))) { next }
    dname_bb <- bs_names$d[bb]
    if (verbose>0) { cat(bs_names$v[bb], "analysis:\n") }

    ## `BayesGLM` call. --------------------------------------------------------
    BayesGLM_results[[dname_bb]] <- BayesGLM(
      BOLD = BOLD[[dname_bb]],
      design = design,
      nuisance = nuisance,
      spatial = spatial[[dname_bb]],
      scale_BOLD = scale_BOLD,
      Bayes = do$Bayesian,
      #EM = do_EM,
      ar_order = ar_order,
      ar_smooth = ar_smooth,
      aic = aic,
      n_threads = n_threads,
      return_INLA = return_INLA,
      verbose = verbose,
      meanTol = meanTol,
      varTol = varTol#,
      #emTol=emTol
    )
  }

  ## Construct beta estimates as `xifti` objects. ------------------------------
  if (verbose>0) cat("Formatting results.\n")

  NULL

  # if (design_type != "multi") {
    estimate_xii <- list(Bayes=NULL, classical=NULL)
    for (method in c("classical", "Bayes")[seq(1+do$Bayesian)]) {
      estimate_xii[[method]] <- BayesGLM_cifti_format_results(
        BayesGLM_results = BayesGLM_results,
        do = do,
        spatial = spatial,
        submeta = submeta,
        session_names = session_names,
        field_names = field_names,
        method = method
      )
    }
  #   bestmodel_xii <- sigma2_xii <- NULL

  # } else {
  #   x <- BayesGLM_cifti_format_results_multi(
  #     BayesGLM_results = BayesGLM_results,
  #     session_names = session_names
  #   )
  #   estimate_xii <- x$estimate_xii
  #   bestmodel_xii <- x$bestmodel_xii
  #   sigma2_xii <- x$sigma2_xii
  #   rm(x)

  #   # [TO DO] get rid
  #   # #stuff we don't have when fitting multiple models
  #   # HRFs <- FIR <- design_FIR <- stimulus <- NULL
  #   # field_names <- field_names
  # }

  result <- list(
    estimate_xii = estimate_xii,
    #bestmodel_xii = bestmodel_xii,
    #sigma2_xii = sigma2_xii,
    nuisance = nuisance,
    #FIR = FIR,
    #design_FIR = design_FIR,
    field_names = field_names,
    session_names = session_names,
    dim = c(n_sess = nS, n_time = nT, n_location = nV0),
    BayesGLM_results = BayesGLM_results
  )
  class(result) <- "BayesGLM_cifti"

  result
}
