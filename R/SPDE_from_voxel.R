#' SPDE from voxel model
#'
#' @param spatial See \code{BayesGLM0}.
#' @return List
#' @keywords internal
SPDE_from_voxel <- function(spatial){
  labels <- spatial$labels
  mask <- labels!=0
  nbhd_order <- spatial$nbhd_order
  buffer <- spatial$buffer

  #[TO DO] Allow the user to additionally specify a mask input excluding certain within-ROI locations
  #Simply remove those locations (in addition to bad data locations, as below) from the labels array before creating SPDE

  ROIs <- unique(labels[mask])
  nR <- length(ROIs)

  #construct the C and G for the SPDE by block-diagonalizing over ROIs
  C_list <- G_list <- spde_list <- vector('list', length=nR)
  for (rr in seq(nR)) {
    mask_rr <- (labels == ROIs[rr])
    # [stop] this breaks [TO DO]
    spde_list[[rr]] <- vol2spde(mask_rr, nbhd_order=nbhd_order, buffer=buffer)
    C_list[[rr]] <- spde_list[[rr]]$mats$C
    G_list[[rr]] <- spde_list[[rr]]$mats$G
  }
  C_sub <- Matrix::bdiag(C_list)
  G_sub <- Matrix::bdiag(G_list)

  #construct the SPDE
  Elog.kappa <- Elog.tau <- 0 #prior means for log(kappa) and log(tau)
  Qlog.kappa <- Qlog.tau <- 0.1 #prior precisions for log(kappa) and log(tau)
  spde <- INLA::inla.spde2.generic(
    M0 = C_sub,
    M1 = G_sub,
    M2 = G_sub%*%solve(C_sub, G_sub),
    theta.mu = c(Elog.kappa, Elog.tau),
    theta.Q = diag(c(Qlog.kappa, Qlog.tau)),
    B0 = matrix(c(0, 1, 0), 1, 3),
    B1 = 2*matrix(c(0, 0, 1), 1, 3),
    B2 = 1
  )

  # [TO DO] test this code for multiple regions, it might break
  # Get indices of data locations.
  data_loc_list <- lapply(spde_list, function(x) which(x$idx2 %in% x$idx))
  data_loc <- data_loc_list[[1]]
  if (nR > 1) {
    before <- 0 # Cumulative sum of previous regions.
    for (r in seq(2,nR)) {
      before <- 0 + nrow(spde_list[[r-1]]$mats$C)
      data_loc_r <- data_loc_list[[r]]
      data_loc_r <- data_loc_r + before
      data_loc <- c(data_loc, data_loc_r)
    }
  }

  # Add `buffer_mask` to `spatial`.
  spatial$buffer_mask <- seq(spde$n.spde) %in% data_loc

  list(
    spde = spde,
    spatial = spatial
  )
}
