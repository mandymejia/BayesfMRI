#' multiGLM for CIFTI
#'
#' Performs classical Bayesian GLM for task fMRI activation with CIFTI-format
#'  data, evaluating multiple design matrices. Includes the pre-processing
#'  steps of nuisance regression. Supports single-session analysis only.
#'
#' @inheritSection Connectome_Workbench_Description Connectome Workbench Requirement
#'
#' @inheritParams BOLD_Param_BayesGLM
#' @param design A 3D numeric array that is locations by fields by designs.
#' @inheritParams TR_Param_BayesGLM
#' @inheritParams brainstructures_Param_BayesGLM
#' @inheritParams resamp_res_Param_BayesGLM
#' @inheritParams nuisance_Param_BayesGLM
#' @inheritParams hpf_Param_BayesGLM
#' @param design_canonical TO DO
#' @inheritParams verbose_Param
#' @inheritParams mean_var_Tol_Param
#'
#' @return An object of class \code{"mGLM"}: a list with elements
#'  \describe{
#'    \item{brainstructures}{\code{data.frame} summarizing the spatial features of each brain structure modeled.}
#'    \item{fields}{\code{data.frame} with the \code{name}, related \code{task}, and \code{HRF_order} of each field.}
#'  }
#'
# @importFrom ciftiTools read_cifti resample_gifti as.xifti remove_xifti
#' @import ciftiTools
#' @importFrom fMRItools match_input is_1
#'
#' @export
#'
multiGLM <- function(
  BOLD,
  design,
  brainstructures=c('left','right'),
  TR=NULL,
  # For surface models
  resamp_res=10000,
  # Nuisance
  hpf=NULL,
  # Below arguments shared with `multiGLM`.
  nuisance=NULL,
  design_canonical=NULL,
  verbose = 1,
  meanTol = 1e-6,
  varTol = 1e-6#,
  #snrTol = 50
){

  scale_design <- FALSE

  # Argument checks. -----------------------------------------------------------
  ### Simple parameters. -------------------------------------------------------
  stopifnot(is.null(TR) || (fMRItools::is_1(TR, "numeric") && TR>0))

  do <- vector("list")

  # In a separate function because these checks are shared with `fit_bayesglm`.
  BayesGLM_argChecks(
    Bayes=FALSE,
    verbose = verbose,
    meanTol = meanTol,
    varTol = varTol
  )

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
  if (nS > 1) { stop("Not supported: multi-session analysis in `multiGLM`.") }

  ### Check `design`. ----------------------------------------------------------
  # Make `design` a sessions-length list of design matrices.
  #   Get `nK`, `field_names`, and `do$perLocDesign`. Check for consistent dims
  #   across sessions.
  x <- BayesGLM_format_design(design, scale_design=FALSE, nS_expect=nS)
  design <- x$design
  nT <- x$nT
  nK <- x$nK
  nD <- x$nD
  field_names <- x$field_names
  design_names <- x$design_names
  do$perLocDesign <- x$per_location_design
  stopifnot(!do$perLocDesign)
  rm(x)

  if (verbose>0) {
    cat("Number of timepoints:    ",
      if (length(unique(nT))==1) { nT[1] } else { paste0(min(nT), "-", max(nT)) }, "\n")
    cat("Number of fields:        ", nK, "\n")
    cat("Brain structures:        ", paste0(brainstructures, collapse=", "), "\n")
    cat("Field names:             ", paste0(field_names, collapse=", "), "\n")
  }

  ### Check `nuisance`. --------------------------------------------------------
  if (!is.null(nuisance)) {
    nuisance <- BayesGLM_format_nuisance(nuisance, nS_expect=nS, nT_expect=nT)
  } else {
    nuisance <- list(NULL)
  }

  ### Make DCT bases in `design` for the high-pass filter. ---------------------
  if (!is.null(hpf)) {
    stopifnot(fMRItools::is_1(hpf, "numeric") && hpf>0)
    DCTs <- lapply(nT, function(nT_ss){
      fMRItools::dct_bases(nT_ss, round(dct_convert(T_=nT_ss, TR=TR, f=hpf)))
    })
    nDCTs <- vapply(DCTs, ncol, 0)
    if (verbose > 0) {
      cat("Including",
        if (length(unique(nDCTs))==1) { nDCTs[1] } else { cat(min(nDCTs), "-", max(nDCTs)) },
        "DCT bases in `nuisance` for highpass filtering.\n")
    }
    nuisance[[1]] <- cbind2(nuisance[[1]], DCTs[[1]] )
  }

  # `BOLD`: read in data and metadata, for each brain structure. ---------------
  nV_T <- setNames(NA*vector("numeric", length(brainstructures)), names(brainstructures))
  names(nV_T)[names(nV_T)=="cortex_left"] <- "cortexL"
  names(nV_T)[names(nV_T)=="cortex_right"] <- "cortexR"

  BOLD_input_msg <- function(do=c("read", "resample")){
    do <- switch(do, read="Reading", resample="Resampling")
    out <- paste0("\t", do,  " BOLD data")
    if (do=="resample") { out <- paste0(out, " to ", resamp_res) }
    out <- paste0(out, ".\n")
  }

  # Above code based on `BayesGLM` which allows `nS>1`. Simplify now.
  BOLD <- BOLD[[1]]
  design <- design[[1]]
  nuisance <- nuisance[[1]]
  ss <- 1

  ### Read and/or resample the CIFTI data. -----------------------------------
  if (is_xifti) {
    if (do$cortex) {
      if (!is.null(resamp_res)) {
        if (any(ciftiTools::infer_resolution(BOLD)!=resamp_res)) {
          if (verbose>0) { cat(BOLD_input_msg("resample")) }
          BOLD <- resample_xifti(BOLD, resamp_res=resamp_res, verbose=FALSE)
        }
      }
    }
  } else {
    if (verbose>0) { cat(BOLD_input_msg("read")) }
    BOLD <- read_cifti(
      BOLD, brainstructures=brainstructures,
      resamp_res=resamp_res
    )
  }

  if (ncol(BOLD) != nT) { stop(
    "The design indicates ", nT,
    " volumes, but the `xifti` data has ", ncol(BOLD),
    " volumes. These must match. Correct either `design` or `BOLD`."
  )}

  if (do$left && is.null(BOLD$data$cortex_left)) { stop("Left cortex data is missing from this BOLD session.") }
  if (do$right && is.null(BOLD$data$cortex_right)) { stop("Right cortex data is missing from this BOLD session.") }
  if (do$sub && is.null(BOLD$data$subcort)) { stop("Subcortex data is missing from this BOLD session.") }

  ### Check `BOLD` data dimensions. ------------------------------------------
  # `xii_res`: total spatial dims according to the BOLD data.
  #   Names: 'left', 'right', 'subcort'.
  xii_res <- vector("numeric", 0)
  if (do$cortex) {
    xii_res <- c(xii_res, ciftiTools::infer_resolution(BOLD))
  }
  if (do$sub) {
    xii_res <- c(xii_res, c(subcort=sum(BOLD$meta$subcort$mask)))
  }

  # Set `nV_T` based on the `xii_res`.
  if (do$left) { nV_T["cortexL"] <- xii_res["left"] }
  if (do$right) { nV_T["cortexR"] <- xii_res["right"] }
  if (do$sub) { nV_T["subcort"] <- sum(BOLD$meta$subcort$mask) }

  rm(xii_res)

   ### Collect metadata. ---------------------------------------------
   # Cortex: ROI `maskL` or `maskR`.
   # Subcortex: `submeta`.`
  if (do$left) {
    maskL <- BOLD$meta$cortex$medial_wall_mask$left
    if (is.null(maskL)) { maskL <- rep(TRUE, nV_T["cortexL"]) }
  } else {
    maskL <- NULL
  }

  if (do$right) {
    maskR <- BOLD$meta$cortex$medial_wall_mask$right
    if (is.null(maskR)) { maskR <- rep(TRUE, nV_T["cortexR"]) }
  } else {
    maskR <- NULL
  }

  if (do$sub) {
    submeta <- BOLD$subcort
  } else {
    submeta <- NULL
  }

  # Collate `BOLD` by brainstructure.
  BOLD <- BOLD$data
  names(BOLD) <- c("cortexL", "cortexR", "subcort")
  if (!do$left) { BOLD$cortexL <- NULL }
  if (!do$right) { BOLD$cortexR <- NULL }
  if (!do$sub) { BOLD$subcort <- NULL }

  # Do GLM. --------------------------------------------------------------------
  mGLM0s <- setNames(vector("list", length(BOLD)), names(BOLD))

  ## Loop through brainstructures. ---------------------------------------------
  bs_names <- data.frame(
    d = c("cortexL", "cortexR", "subcort"), # Names used in this code.
    v = c("Left cortex", "Right cortex", "Subcortex") # Verbose names to show user.
  )

  for (bb in seq(nrow(bs_names))) {
    if (!(bs_names$d[bb] %in% names(BOLD))) { next }
    dname_bb <- bs_names$d[bb]
    if (verbose>0) { cat(bs_names$v[bb], "analysis:\n") }

    ## `multiGLM_fun` call. ----------------------------------------------------
    mGLM0s[[dname_bb]] <- multiGLM_fun(
      BOLD = BOLD[[dname_bb]],
      design = design,
      nuisance = nuisance,
      design_canonical = design_canonical,
      verbose = verbose,
      meanTol = meanTol,
      varTol = varTol
    )
  }

  ## Construct beta estimates as `xifti` objects. ------------------------------
  if (verbose>0) cat("Formatting results.\n")

  result <- list(
    bestmodel_xii = NULL,
    pvalF_xii = NULL,
    Fstat_xii = NULL,
    nuisance = nuisance,
    field_names = field_names,
    mGLM0s = mGLM0s
  )

  #format as xii: (1) index of best model and (2) locations of no signal
  for(meas in c('bestmodel', 'Fstat', 'pvalF')){
    xii_meas <- as.xifti(
      cortexL = mGLM0s$cortexL[[meas]],
      cortexL_mwall = maskL,
      cortexR = mGLM0s$cortexR[[meas]],
      cortexR_mwall = maskR,
      subcortVol = mGLM0s$subcortical[[meas]],
      subcortLabs = submeta$labels,
      subcortMask = submeta$mask
    )
    xii_meas$meta$subcort$trans_mat <- submeta$trans_mat
    xii_meas$meta$subcort$trans_units <- submeta$trans_units
    xii_meas$meta$cifti$names <- field_names
    result[[paste0(meas, "_xii")]] <- xii_meas
  }

  # Convert `bestmodel_xii` to `dlabel`.
  result$bestmodel_xii$meta$cifti$names <- "multiGLM"
  result$bestmodel_xii <- convert_xifti(
    result$bestmodel_xii, "dlabel", 
    levels=seq(nD), levels_old=seq(nD),
    labels=design_names, add_white=FALSE
  )

  class(result) <- "mGLM"
  result
}
