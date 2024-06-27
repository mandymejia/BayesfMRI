#' BayesGLM for CIFTI
#'
#' Performs spatial Bayesian GLM for task fMRI activation with CIFTI-format
#'  data. The cortex is modeled as a surface mesh, and subcortical structures
#'  are modeled as distinct volumetric regions. Includes the pre-processing
#'  steps of nuisance regression, prewhitening, scaling, and variance
#'  normalization. Supports both single- and multi-session analysis. Can also
#'  compute just the classical (spatially-independent)
#'
#' To use \code{BayesGLM}, the design matrix must first be constructed
#'  with \code{\link{make_design}}.
#'
#' @inheritSection Connectome_Workbench_Description Connectome Workbench Requirement
#' @inheritSection INLA_Description INLA Requirement
#' @inheritSection INLA_Latent_Fields_Limit_Description INLA Latent Fields Limit
#'
#' @inheritParams BOLD_Param_BayesGLM
#' @inheritParams brainstructures_Param_BayesGLM
#' @param subROI Which subcortical ROIs should be analyzed?
#' @inheritParams design_Param_BayesGLM
#' @inheritParams nuisance_Param_BayesGLM
#' @inheritParams hpf_Param_BayesGLM
#' @inheritParams TR_Param_BayesGLM
#' @inheritParams surfaces_Param_BayesGLM
#' @inheritParams resamp_res_Param_BayesGLM
#' @inheritParams nbhd_order_Param
#' @inheritParams buffer_Param
#' @inheritParams session_names_Param
#' @inheritParams scale_BOLD_Param
#' @inheritParams Bayes_Param
#' @param hyperpriors Should informative or default non-informative hyperpriors be assumed on SPDE hyperparameters?
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
#' @return An object of class \code{"BayesGLM"}: a list with elements
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
#' @importFrom fMRItools match_input is_1
#' @importFrom matrixStats rowVars rowSums2 colVars
#' @importFrom parallel detectCores
#' @importFrom Matrix bdiag
#'
#' @export
#'
BayesGLM <- function(
  BOLD,
  brainstructures=c("left", "right"),
  subROI=c('Amygdala-L','Amygdala-R','Caudate-L','Caudate-R','Hippocampus-L','Hippocampus-R','Thalamus-L','Thalamus-R'),
  design,
  # Nuisance
  nuisance=NULL,
  hpf=NULL,
  TR=NULL,
  # For surface models
  surfL=NULL,
  surfR=NULL,
  resamp_res=10000,
  # For volume model
  nbhd_order=1,
  buffer=c(1,1,3,4,4),
  # Below arguments shared with `bayesglm_fun`.
  session_names=NULL,
  scale_BOLD = c("mean", "sd", "none"),
  Bayes = TRUE,
  hyperpriors = c("informative","default"),
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
  scale_design <- FALSE # this function no longer does any design matrix construction/alteration besides centering

  # Argument checks. -----------------------------------------------------------
  ### Simple parameters. -------------------------------------------------------
  stopifnot(is.null(TR) || (fMRItools::is_1(TR, "numeric") && TR>0))
  if(Bayes==FALSE) nbhd_order <- 0 #no need for boundary layers with classical GLM
  stopifnot(is.numeric(nbhd_order))
  stopifnot(fMRItools::is_1(nbhd_order, "numeric"))
  stopifnot(nbhd_order>=0 && nbhd_order==round(nbhd_order))
  stopifnot(is.numeric(buffer) || is.null(buffer))

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

  # `hpf` and `TR`.
  if (is.null(hpf)) {
    warn_msg <- if (is.null(TR)) { "`hpf` and `TR`" } else { "`hpf`" }
    warning("Highpass-filtering (HPF) is recommended for computing a GLM on ",
      "time series data, such as fMRI. Set ", warn_msg, " to enable the HPF. ",
      "Or, set `hpf='already'` if the data, design, and nuisance inputs have ",
      "already been high-pass filtered.")
  } else {
    if (fMRItools::is_1(hpf, "character") && hpf=="already") {
      hpf <- NULL
    } else if (is.null(TR)) {
      stop("`hpf` requires `TR`.")
    }
  }

  # Check `BOLD` w/o reading CIFTIs in; check `design` and `nuisance`. ---------
  #   Get all dimensions except for `nV` (because `BOLD` is not read in yet.)

  ### Check `BOLD`. ------------------------------------------------------------
  # Make `BOLD` a sessions-length character vector, or a sessions-length list of
  #  \code{"xifti"} objects. Get `nS`. Do not read or check dims yet.
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

  nS <- length(BOLD)
  if (verbose>0) {
    cat("Number of BOLD sessions: ", nS, "\n")
  }

  ### Brain structures. --------------------------------------------------------
  if ("all" %in% brainstructures) brainstructures <- c("left","right","subcortical")
  brainstructures <- fMRItools::match_input(
    brainstructures, c("left","right","subcortical"),
    user_value_label="brainstructures"
  )
  if (is_xifti) {
    has_bs <- c("left", "right", "subcortical")[!vapply(BOLD[[1]]$data, is.null, FALSE)]
    if(!all(brainstructures %in% has_bs)) stop("BOLD data does not contain all of the structures indicated in `brainstructures`")
  }
  
  do$left <- ('left' %in% brainstructures)
  do$right <- ('right' %in% brainstructures)
  do$sub <- ('subcortical' %in% brainstructures)
  do$cortex <- do$left || do$right
  if (!do$cortex) { resamp_res <- NULL }

  ### Check `design`. ----------------------------------------------------------
  # Make `design` a sessions-length list of design matrices.
  #   Get `nK`, `field_names`, and `do$perLocDesign`. Check for consistent dims
  #   across sessions.
  x <- BayesGLM_format_design(design, scale_design=FALSE, nS_expect=nS)
  design <- x$design
  nT <- x$nT #number of BOLD time points (can vary by session)
  nK <- x$nK #number of regressors in design matrix (fixed across sessions)
  nD <- x$nD #number of design matrices: either 1 or # of BOLD locations
  field_names <- x$field_names
  design_names <- x$design_names
  do$perLocDesign <- x$per_location_design
  rm(x)

  if (do$Bayesian && nK > 5) {
    message(
      "The number of regressors to be modeled spatially exceeds five. ",
      "INLA computation may be slow. Consider reducing the number of design ",
      "matrix columns, e.g. by moving HRF derivative columns to nuisance."
    )
    Sys.sleep(5)
  }

  ### Get `session_names`. -----------------------------------------------------
  session_names <- BayesGLM_session_names(
    nS, session_names, names(BOLD), names(design)
  )
  names(BOLD) <- session_names
  names(design) <- session_names

  if (verbose>0) {
    cat("Session names:           ", paste0(session_names, collapse=", "), "\n")
    cat("Number of timepoints:    ",
        if (length(unique(nT))==1) { nT[1] } else { paste0(min(nT), "-", max(nT)) }, "\n")
    cat("Number of fields:        ", paste0(nK, " (", paste0(field_names, collapse=", "), ")\n"))
  }

  ### Check `nuisance`. --------------------------------------------------------
  if (!is.null(nuisance)) {
    nuisance <- BayesGLM_format_nuisance(nuisance, nS_expect=nS, nT_expect=nT)

    if (!is.null(names(nuisance)) && !all(names(nuisance) == session_names)) {
      #warning("Ignoring `names(nuisance)`; use `session_names` in `make_design`.")
    }
  } else {
    nuisance <- vector("list", nS)
  }
  names(nuisance) <- session_names
  nK2 <- ncol(nuisance[[1]])
  if(is.null(nK2)) nK2 <- 0

  if (verbose>0) {
    cat("Num. nuisance regressors:", nK2, "\n")
  }


  ### Make DCT bases for the high-pass filter ----------------------------------
  if (!is.null(hpf)) {
    stopifnot(fMRItools::is_1(hpf, "numeric") && hpf>0)
    DCTs <- lapply(nT, function(nT_ss){
      fMRItools::dct_bases(nT_ss, round(fMRItools::dct_convert(T_=nT_ss, TR=TR, f=hpf)))
    })
    nDCTs <- vapply(DCTs, ncol, 0)
    if (verbose > 0) {
      cat("Including",
        if (length(unique(nDCTs))==1) { nDCTs[1] } else { cat(min(nDCTs), "-", max(nDCTs)) },
        "DCT bases in `nuisance` for highpass filtering.\n")
    }
    for (ss in seq(nS)) {
      colnames(DCTs[[ss]]) <- paste0('DCT', 1:ncol(DCTs[[ss]]))
      nuisance[[ss]] <- cbind2(nuisance[[ss]], DCTs[[ss]] ) }
  }

  ### Identify missing tasks in design -----------------------------------------

  design_type <- if (do$perLocDesign) { "per_location" } else { "regular" }
  valid_cols <- do.call(rbind, lapply(design, function(q) {
    apply(q, 2, function(r){!all(is.na(r))})
  }))
  if (any(colSums(valid_cols)==0)) { stop("Some tasks are missing (NA) from every session. Please fix.") }
  for (ss in seq(nS)) {
    vcols_ss <- valid_cols[ss,]
    any_bad_design_cols <- if (nD == 1) {
      any(is.na(c(design[[ss]][,vcols_ss])))
    } else {
      any(is.na(c(design[[ss]][,vcols_ss,])))
    }
    if (any_bad_design_cols) {
      stop("`design` has some sessions & tasks for which some data values ",
        "are `NA`. Partially missing data is not allowed. (Missing tasks ",
        "should have all `NA`.)")
    }
  }

  ### Check for intercept (must only be in nuisance) ---------------------------

  des_has_intercept <- vector("logical", nS)
  for (ss in seq(nS)) {
    # Stop if any zero-var, zero-mean column exists.
    des_ss_is_flat <- apply(abs(design[[ss]]) < 1e-8, 2, all)
    if (any(des_ss_is_flat)) {
      stop("Design matrix for session ", ss, " has at least one column that is ",
           "flat (all values are near-zero).")
    }

    # Detect zero-var, nonzero-mean columns.
    if (design_type=="per_location") {
      des_ss_is_intercept <- apply(design[[ss]], c(2,3), var) < 1e-8
    } else {
      des_ss_is_intercept <- matrixStats::colVars(design[[ss]]) < 1e-8
    }

    if(any(des_ss_is_intercept)) {
      stop('Design matrix must not have an intercept, since data and design will ',
           'be centered prior to model fitting. Please fix.')
    }

    # For nuisance: detect zero-var, nonzero-mean columns.
    if (!is.null(nuisance[[ss]])) {
      nuis_ss_is_intercept <- matrixStats::colVars(nuisance[[ss]]) < 1e-8
    } else {
      nuis_ss_is_intercept <- FALSE
    }

    # If no intercept in nuisance, add one.
    if (!any(nuis_ss_is_intercept)) {
      nuisance[[ss]] <- cbind2(nuisance[[ss]], as.matrix(rep(1, nT[ss])))
    }

    ### CHECK FOR COLLINEARITY AMONG ALL REGRESSORS (DESIGN + NUISANCE + DCT)

    #Check for collinearity among columns of design and nuisance
    nuis_ss <- nuisance[[ss]]
    checkCorr <- function(design){
      X <- cbind(nuis_ss, design)
      suppressWarnings(corX <- cor(X)); diag(corX) <- NA
      corX
    }
    checkVIF <- function(design){
      int <- (apply(nuis_ss, 2, var) == 0) #exclude intercept column of nuisance
      X <- cbind(nuis_ss[,!int], design)
      y <- rnorm(nT[ss]) #add fake y variable, has no influence
      Xnames <- paste0("X",1:ncol(X))
      df <- as.data.frame(cbind(X, y)); names(df) <- c(Xnames,"y")
      f <- as.formula(paste0('y ~ ',paste(Xnames, collapse = " + ")))
      suppressWarnings(car::vif(lm(f, data = df)))
    }


    # Single Design Matrix
    vcols_ss <- valid_cols[ss,]
    if (!do$perLocDesign) {
      cor_x <- checkCorr(design[[ss]][,vcols_ss])
      x1 <- max(abs(cor_x), na.rm=TRUE) #max correlation between any two columns of design & nuisance
      x2 <- checkVIF(design[[ss]][,vcols_ss])
      x2a <- x2[1:(ncol(nuis_ss) - 1)] #nuisance regressors
      x2b <- x2[ncol(nuis_ss):length(x2)] #design regressors
      if(verbose > 0) {
        cat('Checking for collinearity of the design matrix and nuisance matrix (including DCT bases) collectively \n')
        cat(paste0('\tVIF for design regressors: ', paste0(round(x2b), collapse=', '),'\n'))
        cat(paste0('\tMaximum VIF among all nuisance regressors: ', round(max(x2a)),'\n'))
        inds <- which(abs(cor_x) == x1, arr.ind = TRUE)
        cat(paste0('\tMaximum correlation among all regressors: ', round(x1,4), ' (',
            rownames(cor_x)[inds[1,1]], ' and ', rownames(cor_x)[inds[1,2]], ')\n'))
      }
#
#       if(verbose > 0) {
#         if(x1 > 0.95) {
#           inds <- which(abs(cor_x) > 0.95, arr.ind = TRUE)
#           rows <- unique(inds[,1])
#           cols <- unique(inds[,2])
#           print(round(cor_x[rows,cols],2))
#         }
#       }

      if(x1 > 0.999) {
        stop('I detected high collinearity (cor > 0.999) between regressors in the design and nuisance matrices. Please fix.')
      }
    } else {
      # Multiple Design Matrices (one per location)
      x1 <- apply(
        design[[ss]][,vcols_ss,,drop=FALSE], 3, function(x)
        max(abs(checkCorr(x)), na.rm=TRUE)
      )
      if(verbose > 0) {
        cat('Checking for collinearity of the design matrix and nuisance matrix (including DCT bases) collectively \n')
        cat(paste0('\tMaximum correlation among regressors, max over locations: ', round(max(x1),2),'\n'))
      }
      if(max(x1) > 0.99) stop('I detected high collinearity (cor > 0.99) between regressors in the design and nuisance matrices for at least one location. Please fix.')
    }
  }

  # Initialize `spatial` to store all spatial information. ---------------------
  spatial <- list(
    cortexL = list(spatial_type="surf", surf=NULL, mask=NULL),
    cortexR = list(spatial_type="surf", surf=NULL, mask=NULL),
    subcort = list(
      spatial_type="voxel",
      labels=NULL,
      trans_mat=NULL, trans_units=NULL,
      nbhd_order=nbhd_order, buffer=buffer,
      buffer_mask=NULL, # created in `SPDE_from_voxel`
      data_loc=NULL # created in `fit_bayesglm`
    )
  )
  if (!do$left) { spatial$cortexL <- NULL }
  if (!do$right) { spatial$cortexR <- NULL }
  if (!do$sub) { spatial$subcort <- NULL }

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
        if(verbose>0) cat("Since no surfaces provided or present, I will use the fs_LR inflated surfaces for all modeling.\n")
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
        if(verbose>0) cat("Since no surfaces provided or present, I will use the fs_LR inflated surfaces for all modeling.\n")
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

  # Total number of locations, `nV_T`. Check it's consistent across sessions.
  #   Just use to check for the same resolution across multiple-sessions.
  #   Instead of passing `nV_T`, `spatial` is passed on to `BayesGLM`.
  nV_T <- setNames(NA*vector("numeric", length(spatial)), names(spatial))

  BOLD_input_msg <- function(ss, nS, do=c("read", "resample")){
    do <- switch(do, read="Reading", resample="Resampling")
    out <- if (ss==1) { paste0("\t", do,  " BOLD data") } else { "" }
    if (do=="resample") { out <- paste0(out, " to ", resamp_res) }
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
            BOLD[[ss]] <- resample_xifti(BOLD[[ss]], resamp_res=resamp_res, verbose=FALSE)
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

    #Check that all brainstructures are present
    has_bs <- c("left", "right", "subcortical")[!vapply(BOLD[[ss]]$data, is.null, FALSE)]
    if(!all(brainstructures %in% has_bs)){
      if(nS == 1) stop("BOLD data does not contain all of the structures indicated in `brainstructures`")
      if(nS > 1) stop(paste0("Session ", ss, " BOLD data does not contain all of the structures indicated in `brainstructures`"))
    }

    #Remove extra subcortical ROIs
    if(do$sub){
      label_names <- as.character(unique(BOLD[[ss]]$meta$subcort$labels))
      if(!all(subROI %in% label_names)) stop('All elements of subROI must be valid subcortical labels and present in BOLD')
      mask_new <- BOLD[[ss]]$meta$subcort$labels %in% subROI
      BOLD[[ss]]$data$subcort <- BOLD[[ss]]$data$subcort[mask_new,]
      BOLD[[ss]]$meta$subcort$mask[BOLD[[ss]]$meta$subcort$mask] <- mask_new
      BOLD[[ss]]$meta$subcort$labels <- BOLD[[ss]]$meta$subcort$labels[mask_new]
    }

    if (ss == 1 && verbose>0) {
      cat("Brain structures:   ", paste0(brainstructures, collapse=", "), "\n")
      if (do$sub) {
        cat("Subcortical ROIs: ", paste0(subROI, collapse=", "), "\n")
      }
    }

    #Check that nT matches design matrix
    if (ncol(BOLD[[ss]]) != nT[ss]) { stop(
      "The design for session '", session_names[ss], "' indicates ", nT[ss],
      " volumes, but the `xifti` data for this session has ", ncol(BOLD[[ss]]),
      " volumes. These must match. Correct either `design` or `BOLD`."
    )}

    if (do$perLocDesign) {
      if (nrow(BOLD[[ss]]) != nD) { stop(
        "The design matrix array indicates ", nD,
        " total locations, each being modeled with its own design matrix. ",
        "However, the `xifti` data for this session has ", nrow(BOLD[[ss]]),
        " total locations. These must match. Correct either `design` or `BOLD`."
      )}
    }

    #should be redundant with above
    if (do$left && is.null(BOLD[[ss]]$data$cortex_left)) { stop("Left cortex data is missing from this BOLD session.") }
    if (do$right && is.null(BOLD[[ss]]$data$cortex_right)) { stop("Right cortex data is missing from this BOLD session.") }
    if (do$sub && is.null(BOLD[[ss]]$data$subcort)) { stop("Subcortex data is missing from this BOLD session.") }

    ### Check `BOLD` data dimensions. ------------------------------------------
    # `xii_res`: total spatial dims according to the BOLD data.
    #   Names: 'left', 'right', 'subcort'.
    xii_res <- vector("numeric", 0)
    if (do$cortex) {
      xii_res <- c(xii_res, ciftiTools::infer_resolution(BOLD[[ss]]))
    }
    if (do$sub) {
      xii_res <- c(xii_res, c(subcort=sum(BOLD[[ss]]$meta$subcort$mask)))
    }

    # Set `nV_T` based on the `xii_res` of the first session...
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
        nV_T["cortexL"] <- xii_res["left"]
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
        nV_T["cortexR"] <- xii_res["right"]
      }
      if (do$sub) {
        nV_T["subcort"] <- sum(BOLD[[ss]]$meta$subcort$mask)
      }

    # ...Check `nV_T` matches `xii_res` of other sessions.
    } else {
      if (do$left && (xii_res["left"] != nV_T["cortexL"])) {
        stop(
          "This BOLD session appears to have left cortex data with ",
          xii_res["left"], "total vertices, while the first session ",
          "appears to have", nV_T["cortexL"], "total vertices."
        )
      }
      if (do$right && (xii_res["right"] != nV_T["cortexR"])) {
        stop(
          "This BOLD session appears to have right cortex data with ",
          xii_res["right"], "total vertices, while the first session ",
          "appears to have", nV_T["cortexR"], "total vertices."
        )
      }
      if (do$sub && (xii_res["subcort"] != nV_T["subcort"])) {
        stop(
          "This BOLD session appears to have subcortical data with",
          xii_res["subcort"], "total voxels, while the first session ",
          "appears to have", nV_T["subcort"], " total voxels."
        )
      }
    }
    rm(xii_res)
  }

   ### Collect `spatial` metadata. ---------------------------------------------
   # Cortex: ROI `mask`. Use the intersection.
   # Subcortex: `labels`, `trans_mat`, `trans_units`. Require same across sessions.
  if (do$left) {
    maskL <- lapply(BOLD, function(q){
      mask_ss <- q$meta$cortex$medial_wall_mask$left
      if (is.null(mask_ss)) { mask_ss <- rep(TRUE, nV_T["cortexL"]) }
      mask_ss
    })
    maskL <- colSums(do.call(rbind, maskL))
    maskL_has_diffs <- (!all(maskL %in% c(0, nS)))
    if (maskL_has_diffs && verbose>0) {
      cat("BOLD left cortex initial ROIs do not match across sessions; using their intersection by setting mismatch locations to `NA`.\n")
    }
    spatial$cortexL$mask <- maskL == nS
    if (maskL_has_diffs) {
      # Set to `NA`--masking of data (and design for per-location model) will be handled in `BayesGLM`.
      for (ss in seq(nS)) {
        new_mask_ss <- maskL[BOLD[[ss]]$meta$cortex$medial_wall_mask$left]
        BOLD[[ss]]$data$cortex_left[!new_mask_ss,] <- NA
      }
    }
    #if(verbose > 0) cat('Left cortex resolution: ', paste0(nrow(BOLD[[1]]$data$cortex_left),' vertices\n'))
    rm(maskL, maskL_has_diffs)
  }

  if (do$right) {
    maskR <- lapply(BOLD, function(q){
      mask_ss <- q$meta$cortex$medial_wall_mask$right
      if (is.null(mask_ss)) { mask_ss <- rep(TRUE, nV_T["cortexR"]) }
      mask_ss
    })
    maskR <- colSums(do.call(rbind, maskR))
    maskR_has_diffs <- (!all(maskR %in% c(0, nS)))
    if (maskR_has_diffs && verbose>0) {
      cat("BOLD right cortex initial ROIs do not match across sessions; using their intersection by setting mismatch locations to `NA`.\n")
    }
    spatial$cortexR$mask <- maskR == nS
    if (maskR_has_diffs) {
      # Set to `NA`--masking of data (and design for per-location model) will be handled in `BayesGLM`.
      for (ss in seq(nS)) {
        new_mask_ss <- maskR[BOLD[[ss]]$meta$cortex$medial_wall_mask$right]
        BOLD[[ss]]$data$cortex_right[!new_mask_ss,] <- NA
      }
    }
    #if(verbose > 0) cat('Right cortex resolution: ', paste0(nrow(BOLD[[1]]$data$cortex_right),' vertices\n'))
    rm(maskR, maskR_has_diffs)
  }


  if (do$sub) {
    for (ss in seq(nS)) {
      if (ss == 1) {
        submeta <- BOLD[[ss]]$meta$subcort
        spatial$subcort["labels"] <- list(submeta$mask*0)
        spatial$subcort$labels[submeta$mask==TRUE] <- submeta$labels
        spatial$subcort["trans_mat"] <- list(submeta$trans_mat)
        spatial$subcort["trans_units"] <- list(submeta$trans_units)
      } else {
        stopifnot(length(dim(spatial$subcort$labels)) == length(dim(BOLD[[ss]]$meta$subcort$mask)))
        stopifnot(all(dim(spatial$subcort$labels) == dim(BOLD[[ss]]$meta$subcort$mask)))
        stopifnot(all((spatial$subcort$labels!=0) == BOLD[[ss]]$meta$subcort$mask))
      }
    }

    #check and report trans_mat
    res <- abs(diag(spatial$subcort$trans_mat)[1:3])
    if(!is.numeric(res)) stop('I cannot infer subcortical voxel resolution from CIFTI header.  Check trans_mat or contact developer.')
    if(any(is.na(res)) | any(is.nan(res))) stop('I cannot infer subcortical voxel resolution from CIFTI header.  Check trans_mat or contact developer.')
    if(min(res) < 1 | max(res > 4)) stop('Voxel resolution appears to be implausible (less than 1 or greater than 4).  Check trans_mat in CIFTI header or contact developer.')
    if(verbose > 0) cat('Subcortical voxel size: ', paste0(paste(res, collapse = ' x '),' mm \n'))

  } else {
    submeta <- NULL
  }

  # Collate `BOLD` by brainstructure.
  BOLD <- lapply(BOLD, '[[', "data")
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

  # Do GLM. --------------------------------------------------------------------

  BGLMs <- setNames(vector("list", length(spatial)), names(spatial))

  ## Loop through brainstructures. ---------------------------------------------
  bs_names <- data.frame(
    d = c("cortexL", "cortexR", "subcort"), # Names used in this code.
    v = c("Left cortex", "Right cortex", "Subcortex") # Verbose names to show user.
  )

  # number of data locations (vs. `nV_T` includes masked locations on the mesh.)

  nV_D <- vapply(lapply(BOLD, function(q){q[[1]]}), ncol, 0)

  for (bb in seq(nrow(bs_names))) {
    if (!(bs_names$d[bb] %in% names(BOLD))) { next }
    dname_bb <- bs_names$d[bb]
    if (verbose>0) { cat(bs_names$v[bb], "analysis:\n") }
    design_bb <- if (do$perLocDesign) {
      lapply(design, function(q){q[,,seq(
        sum(c(0, nV_D)[seq(bb)])+1, sum(nV_D[seq(bb)])
        ),drop=FALSE]})
    } else {
      design
    }

    ## `fit_bayesglm` call. --------------------------------------------------------
    BGLMs[[dname_bb]] <- fit_bayesglm(
      BOLD = BOLD[[dname_bb]],
      design = design_bb,
      nuisance = nuisance,
      spatial = spatial[[dname_bb]],
      scale_BOLD = scale_BOLD,
      Bayes = do$Bayesian,
      hyperpriors = hyperpriors,
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

  estimate_xii <- RSS_xii <- list(Bayes=NULL, classical=NULL)
  for (method in c("classical", "Bayes")[seq(1+do$Bayes)]) {
    x <- BayesGLM_format_cifti(
      BGLMs = BGLMs,
      do = do,
      spatial = spatial,
      submeta = submeta,
      session_names = session_names,
      field_names = field_names,
      method = method # it's 'Bayesian' here not 'Bayes', but still matches.
    )
    estimate_xii[[method]] <- x$estimates
    RSS_xii[[method]] <- x$RSS
    rm(x)
  }

  # [TO DO] HRF_info

  result_dim <- c(
    c(sess = nS, time = nT),
    setNames(nV_T, paste0("loc_", names(nV_T)))
  )

  result <- list(
    estimate_xii = estimate_xii,
    RSS_xii = RSS_xii,
    nuisance = nuisance,
    field_names = field_names,
    session_names = session_names,
    dim = result_dim,
    BGLMs = BGLMs
  )
  class(result) <- "BGLM"

  result
}
