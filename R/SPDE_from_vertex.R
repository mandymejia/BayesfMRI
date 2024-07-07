#' SPDE from mesh model
#'
#' @param spatial See \code{BayesGLM}.
#' @param qc_mask QC mask.
#' @param logkappa,logtau vector of min, max and initial value for prior on log kappa and log tau. Min and max are extreme quantiles, not hard constraints.
#' @return List: \code{mesh}, \code{spde}, \code{spatial}.
#' @keywords internal
SPDE_from_vertex <- function(spatial, qc_mask, logkappa = NULL, logtau = NULL){
  surf <- spatial$surf
  mask <- spatial$maskIn

  # Create INLA mesh.
  mesh <- make_mesh(surf$vertices, surf$faces)

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
  spde <- INLA::inla.spde2.matern(
    mesh,
    prior.tau = exp(Elog.tau),
    prior.kappa = exp(Elog.kappa),
    theta.prior.prec = diag(c(Qlog.tau, Qlog.kappa))
  )

  spatial$maskMdat <- spatial$maskIn; spatial$maskMdat[spatial$maskMdat][!qc_mask] <- FALSE
  spatial$maskMbuf <- !spatial$maskMdat
  spatial$Mmap <- which(spatial$maskMdat)

  list(
    mesh = mesh,
    spde = spde,
    spatial = spatial
  )
}
