#' SPDE from voxel model
#'
#' @param spatial See \code{BayesGLM}.
#' @param qc_mask QC mask.
#' @param logkappa,logtau vector of min, max and initial value for prior on log kappa and log tau. Min and max are extreme quantiles, not hard constraints.
#' @return List: \code{spde} and \code{spatial}.
#' @keywords internal
SPDE_from_voxel <- function(spatial, qc_mask=NULL, logkappa = NULL, logtau = NULL){

  if (!is.null(qc_mask)) {
    # Update `spatial`.
    ### Update `maskMdat` by applying `qc_mask` to it.
    spatial$maskMdat <- spatial$maskIn
    spatial$maskMdat[spatial$maskMdat][!qc_mask] <- FALSE
    ### Update `labsMdat` by using the new `maskMdat`
    spatial$labsMdat <- spatial$labels[spatial$maskMdat[spatial$maskIn]]
  }

  labels <- spatial$labels
  nbhd_order <- spatial$nbhd_order
  buffer <- spatial$buffer
  res <- abs(diag(spatial$trans_mat)[1:3])

  #[TO DO] Allow the user to additionally specify a mask input excluding certain within-ROI locations
  #Simply remove those locations (in addition to bad data locations) from the labels array before creating SPDE

  ROIs <- unique(labels)
  nR <- length(ROIs)

  #construct the C and G for the SPDE by block-diagonalizing over ROIs
  rr_list <- setNames(vector('list', length=nR), ROIs)
  C_list <- G_list <- spde_list <- mask_In_to_Mdat <- rr_list
  maskIn_to_Mdat_all <- spatial$maskMdat[][spatial$maskIn[]]
  for (rr in seq(nR)) {
    # Ge maskIn for just ROI rr
    mask_rr <- spatial$maskIn
    mask_rr[mask_rr][spatial$labels != ROIs[rr]] <- FALSE
    # Construct the SPDE for ROI rr based on maskIn
    spde_list[[rr]] <- vol2spde(mask_rr, nbhd_order=nbhd_order, buffer=buffer, res=res)
    C_list[[rr]] <- spde_list[[rr]]$mats$C
    G_list[[rr]] <- spde_list[[rr]]$mats$G
    # We will drop the QC-masked locations using `dat_loc`. So for now,
    #   Get the mask to convert from maskIn to maskMdat for just ROI rr.
    mask_In_to_Mdat[[rr]] <- maskIn_to_Mdat_all[spatial$labels == ROIs[rr]]
  }
  C_sub <- Matrix::bdiag(C_list)
  G_sub <- Matrix::bdiag(G_list)

  #construct hyperpriors
  Elog.kappa <- Elog.tau <- 0 #prior means for log(kappa) and log(tau)
  Qlog.kappa <- Qlog.tau <- 0.1 #prior precisions for log(kappa) and log(tau)
  if(!is.null(logtau)){
    Elog.tau <- logtau[3]
    prior_sd <- abs(diff(logtau[1:2]))/4 #so that 95% of the prior density is within the range
    Qlog.tau <- 1/(prior_sd^2)
  }
  if(!is.null(logkappa)){
    Elog.kappa <- logkappa[3]
    prior_sd <- abs(diff(logkappa[1:2]))/4 #so that 95% of the prior density is within the range
    Qlog.kappa <- 1/(prior_sd^2)
  }

  #construct the SPDE
  spde <- INLA::inla.spde2.generic(
    M0 = C_sub,
    M1 = G_sub,
    M2 = G_sub%*%solve(C_sub, G_sub),
    theta.mu = c(Elog.tau, Elog.kappa),
    theta.Q = diag(c(Qlog.tau, Qlog.kappa)),
    B0 = matrix(c(0, 1, 0), 1, 3),
    B1 = 2*matrix(c(0, 0, 1), 1, 3),
    B2 = 1
  )

  spatial$nV_M <- spde$n.spde

  # [TO DO] test this code for multiple regions, it might break
  # Get indices of data locations.
  data_loc_list <- lapply(spde_list, function(x) which(x$idx2 %in% x$idx))
  # Drop QC-masked locations from `data_loc_list`
  for (rr in seq(nR)) {
    data_loc_list[[rr]] <- data_loc_list[[rr]][mask_In_to_Mdat[[rr]]]
  }
  # Get vector of data_loc
  data_loc <- data_loc_list[[1]]
  if (nR > 1) {
    before <- 0 # Cumulative sum of previous regions.
    for (rr in seq(2,nR)) {
      before <- before + nrow(spde_list[[rr-1]]$mats$C)
      data_loc_r <- data_loc_list[[rr]]
      data_loc_r <- data_loc_r + before
      data_loc <- c(data_loc, data_loc_r)
    }
  }

  # Add masks and data locations to `spatial`.
  spatial$Mmap <- data_loc
  # maskM <- seq(spde$n.spde) %in% data_loc

  list(
    spde = spde,
    spatial = spatial
  )
}
