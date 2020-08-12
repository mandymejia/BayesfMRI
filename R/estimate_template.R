#' Estimate Template for Template or Diagnostic ICA
#'
#' @param cifti_fnames Vector of file paths of CIFTI-format fMRI timeseries 
#'  (*.dtseries.nii) for template estimation
#' @param cifti_fnames2 Vector of file paths of "retest" CIFTI-format fMRI 
#'  timeseries (*.dtseries.nii) for template estimation.  Must be from the same 
#'  subjects and in the same order as cifti_fnames.  If none specified, will 
#'  create pseudo test-retest data from single session.
#' @param GICA_fname File path of CIFTI-format group ICA maps (ending in .d*.nii)
#' @param inds Indicators of which group ICs to include in template. If NULL, 
#'  use all group ICs.
#' @param scale Logical indicating whether BOLD data should be scaled by the 
#'  spatial standard deviation before template estimation.
#' @param brainstructures Character vector indicating which brain structure(s) 
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right 
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures). 
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param verbose If TRUE, display progress updates
#'
#' @return
#' @export
#' @importFrom ciftiTools read_cifti
#'
estimate_template.cifti <- function(
  cifti_fnames, cifti_fnames2=NULL, GICA_fname, 
  inds=NULL, scale=TRUE, brainstructures=c('left','right'), verbose=TRUE){

  # Check arguments.
  if (!is.logical(scale) || length(scale) != 1) { stop('scale must be a logical value') }

  brainstructures <- ciftiTools:::match_input(
    brainstructures, c("left","right","subcortical","all"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) { 
    brainstructures <- c("left","right","subcortical")
  }

  # Read GICA result. 
  #'  First, obtain the mapping (used to infer brainstructures)
  #   with `cifti_read_flat` for `cifti_fnames` and `cifti_fnames2`).
  if(verbose) cat('\n Reading in GICA result')
  flat_CIFTI_map <- ciftiTools:::map_cifti(GICA_fname)
  #   Next, read the CIFTI. If cortex data to be included does not have a medial
  #   wall, use `cifti_read_separate` to try to infer it from the NIFTI.
  no_left_mwall <- "left" %in% brainstructures && all(flat_CIFTI_map$cortex$medial_wall_mask$left)
  no_right_mwall <- "right" %in% brainstructures && all(flat_CIFTI_map$cortex$medial_wall_mask$left)
  if (no_left_mwall || no_right_mwall) {
    GICA <- read_cifti(
      GICA_fname,  full_volume=TRUE, brainstructures=brainstructures
    )
    no_left_mwall <- "left" %in% brainstructures && all(GICA$meta$cortex$medial_wall_mask$left)
    no_right_mwall <- "right" %in% brainstructures && all(GICA$meta$cortex$medial_wall_mask$left)
    if (no_left_mwall || no_right_mwall) {
      warning(paste(
        "No medial wall vertices were detected in the",
        c("left cortex", "right cortex", "cortex")[no_left_mwall*1 + no_right_mwall*2],
        "component of the GICA CIFTI."
      ))
    # Replace the empty medial wall mask(s).
    } else {
      if ("left" %in% brainstructures) {
        stopifnot(length(flat_CIFTI_map$cortex$medial_wall_mask$left) == length(GICA$meta$cortex$medial_wall_mask$left))
        flat_CIFTI_map$cortex$medial_wall_mask$left <- GICA$meta$cortex$medial_wall_mask$left
      }
      if ("right" %in% brainstructures) {
        stopifnot(length(flat_CIFTI_map$cortex$medial_wall_mask$right) == length(GICA$meta$cortex$medial_wall_mask$right))
        flat_CIFTI_map$cortex$medial_wall_mask$right <- GICA$meta$cortex$medial_wall_mask$right
      }
    }
  #   Also use the `separate` method if subcortical data is included (for
  #   visualization). Otherwise, use the `convert` method.
  } else {
    GICA <- read_cifti(
      GICA_fname, 
      full_volume="subcortical" %in% brainstructures, 
      brainstructures=brainstructures
    )
  }

  # Obtain the brainstructure mask for the flattened CIFTIs.
  #   It will remove any newly-detected medial wall vertices.
  flat_bs_labs <- c(
    ifelse(flat_CIFTI_map$cortex$medial_wall_mask$left, "left", "mwall"),
    ifelse(flat_CIFTI_map$cortex$medial_wall_mask$right, "right", "mwall"),
    rep("subcortical", length(flat_CIFTI_map$subcort$labels))
  )
  flat_bs_mask <- flat_bs_labs %in% brainstructures

  # Flatten. `GICA_flat` will have the subcortical voxels in alphabetical order
  #   because `cifti_read_flat()` returns them in alphabetical order. 
  GICA_flat <- GICA
  if ("subcortical" %in% brainstructures) {
    alpha_order <- order(GICA_flat$meta$subcort$labels)
    GICA_flat$data$subcort <- GICA_flat$data$subcort[alpha_order,, drop=FALSE]
  }
  GICA_flat <- do.call(rbind, GICA_flat$data)
  V <- nrow(GICA_flat); L0 <- ncol(GICA_flat)
  
  # Center each IC map.
  GICA_flat <- scale(GICA_flat, scale=FALSE) 

  if(verbose){
    cat(paste0('\n Number of data locations: ',V))
    cat(paste0('\n Number of group ICs: ',L0))
  }

  L <- L0
  if(!is.null(inds)){
    if(any(!(inds %in% 1:L0))) stop('Invalid entries in inds argument.')
    L <- length(inds)
  } else {
    inds <- 1:L0
  }

  N <- length(cifti_fnames)

  if(verbose){
    cat(paste0('\n Number of template ICs: ',L))
    cat(paste0('\n Number of training subjects: ',N))
  }

  retest <- !is.null(cifti_fnames2)
  if(retest){
    if(length(cifti_fnames) != length(cifti_fnames2)) stop('If provided, cifti_fnames2 must have same length as cifti_fnames and be in the same subject order.')
  }

  # PERFORM DUAL REGRESSION ON (PSEUDO) TEST-RETEST DATA
  DR1 <- DR2 <- array(NA, dim=c(N, L, V))
  missing_data <- NULL
  for(ii in 1:N){

    if(verbose) cat(paste0('\n Reading in data for subject ',ii,' of ',N))

    #read in BOLD
    fname_ii <- cifti_fnames[ii]
    if(!file.exists(fname_ii)) {
      missing_data <- c(missing_data, fname_ii)
      if(verbose) cat(paste0('\n Data not available'))
      next
    }
    BOLD1_ii <- read_cifti(fname_ii, flat=TRUE)[flat_bs_mask[flat_bs_labs != "mwall"],, drop=FALSE]

    if(nrow(BOLD1_ii) != nrow(GICA2)) stop(paste0('The number of data locations in GICA and timeseries data from subject ',ii,' do not match.'))
    ntime <- ncol(BOLD1_ii)

    #read in BOLD retest data OR create pseudo test-retest data
    if(!retest){
      part1 <- 1:round(ntime/2)
      part2 <- setdiff(1:ntime, part1)
      BOLD2_ii <- BOLD1_ii[,part2]
      BOLD1_ii <- BOLD1_ii[,part1]
    } else {
      #read in BOLD from retest
      fname_ii <- cifti_fnames2[ii]
      if(!file.exists(fname_ii)) {
        missing_data <- c(missing_data, fname_ii)
        if(verbose) cat(paste0('\n Data not available'))
        next
      }
      BOLD2_ii <- read_cifti(fname_ii, flat=TRUE)[flat_bs_mask[flat_bs_labs != "mwall"],, drop=FALSE]
    }

    #perform dual regression on test and retest data
    DR1_ii <- dual_reg(t(BOLD1_ii), t(GICA2), scale=scale)$S
    DR2_ii <- dual_reg(t(BOLD2_ii), t(GICA2), scale=scale)$S
    DR1[ii,,] <- DR1_ii[inds,]
    DR2[ii,,] <- DR2_ii[inds,]
  }

  # ESTIMATE MEAN

  if(verbose) cat('\n Estimating Template Mean')
  mean1 <- apply(DR1, c(2,3), mean, na.rm=TRUE)
  mean2 <- apply(DR2, c(2,3), mean, na.rm=TRUE)
  template_mean <- (mean1 + mean2)/2;

  # ESTIMATE SIGNAL (BETWEEN-SUBJECT) VARIANCE

  # total variance
  if(verbose) cat('\n Estimating Total Variance')
  var_tot1 <- apply(DR1, c(2,3), var, na.rm=TRUE)
  var_tot2 <- apply(DR2, c(2,3), var, na.rm=TRUE)
  var_tot <- (var_tot1 + var_tot2)/2

  # noise (within-subject) variance
  if(verbose) cat('\n Estimating Within-Subject Variance')
  DR_diff = DR1 - DR2;
  var_noise <- (1/2)*apply(DR_diff, c(2,3), var, na.rm=TRUE)

  # signal (between-subject) variance
  if(verbose) cat('\n Estimating Template (Between-Subject) Variance')
  template_var <- var_tot - var_noise
  template_var[template_var < 0] <- 0

  rm(DR1, DR2, mean1, mean2, var_tot1, var_tot2, var_tot, DR_diff)

  # Format template as "xifti"s
  xifti_mean <- GICA; xifti_var <- GICA
  if ("left" %in% brainstructures) {
    xifti_mean$data$cortex_left <- template_mean[flat_bs_labs[flat_bs_mask]=="left",, drop=FALSE]
    xifti_var$data$cortex_left <- template_var[flat_bs_labs[flat_bs_mask]=="left",, drop=FALSE]
  }
  if ("right" %in% brainstructures) {
    xifti_mean$data$cortex_right <- template_mean[flat_bs_labs[flat_bs_mask]=="right",, drop=FALSE]
    xifti_var$data$cortex_right <- template_var[flat_bs_labs[flat_bs_mask]=="right",, drop=FALSE]
  }
  if ("subcortical" %in% brainstructures) {
    xifti_mean$data$subcort <- template_mean[flat_bs_labs[flat_bs_mask]=="subcortical",, drop=FALSE]
    xifti_var$data$subcort <- template_var[flat_bs_labs[flat_bs_mask]=="subcortical",, drop=FALSE]
  }

  list(mean_template=xifti_mean, var_template=xifti_var)
}