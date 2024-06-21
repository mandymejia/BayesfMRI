#' SPDE from mesh model
#'
#' @param spatial See \code{BayesGLM0}.
#' @param logkappa,logtau vector of min, max and initial value for prior on log kappa and log tau. Min and max are extreme quantiles, not hard constraints.
#' @return List
#' @keywords internal
SPDE_from_mesh <- function(spatial, logkappa = NULL, logtau = NULL){
  surf <- spatial$surf
  mask <- spatial$mask
  nV <- nrow(surf$vertices)

  # Create INLA mesh.
  mesh_full <- make_mesh(surf$vertices, surf$faces)

  # [TO DO] change strategy; use a full mesh.
  # Create submesh of in-mask locations.
  # Update masks (sometimes `submesh` will exclude additional vertices).
  if (!all(mask)) {
    mesh <- submesh(mask, mesh_full)
    mask_new <- !is.na(mesh$idx$loc)
    mask_new_diff <- mask_new[mask]
    mask[mask] <- mask_new_diff
    mesh$idx$loc <- mesh$idx$loc[mask]
  }

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

  #construct SPDE
  spde <- INLA::inla.spde2.matern(mesh,
                                  prior.tau = exp(Elog.tau),
                                  prior.kappa = exp(Elog.kappa),
                                  theta.prior.prec = diag(c(Qlog.tau, Qlog.kappa)))

  list(
    mesh_full = mesh_full,
    mesh = mesh,
    mask_new_diff = mask_new_diff,
    spde = spde,
    spatial = list(surf=surf, mask=mask),
    data_loc = NULL # [TO DO] implement this so we can have boundary vertices, i.e. the medial wall
  )
}
