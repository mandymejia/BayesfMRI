#' Calculate the SPDE covariance
#'
#' @param kappa2 A scalar
#' @param phi A scalar
#' @param spde An \code{inla.spde2} object containing the information about the
#'   mesh structure for the SPDE prior
#'
#' @return The SPDE prior matrix
#' @keywords internal
spde_Q_phi <- function(kappa2, phi, spde) {
  Cmat <- spde$param.inla$M0
  Gmat <- (spde$param.inla$M1 + Matrix::t(spde$param.inla$M1))/2
  GtCinvG <- spde$param.inla$M2
  Q <- (kappa2*Cmat + 2*Gmat + GtCinvG/kappa2) / (4*pi*phi)
  return(Q)
}

#' Gives the portion of the Q matrix independent of phi
#'
#' @param kappa2 scalar
#' @param spdean spde object
#'
#' @return a dgCMatrix
#' @keywords internal
Q_prime <- function(kappa2, spde) {
  Cmat <- spde$param.inla$M0
  Gmat <- (spde$param.inla$M1 + Matrix::t(spde$param.inla$M1))/2
  GtCinvG <- spde$param.inla$M2
  Q <- (kappa2*Cmat + 2*Gmat + GtCinvG/kappa2)
  return(Q)
}

#' The negative of the objective function for kappa
#'
#' @param kappa2 scalar
#' @param spde an spde object
#' @param phi scalar
#' @param Sigma dgCMatrix
#' @param mu dgeMatrix
#'
#' @return a scalar
#' @keywords internal
neg_kappa_fn <- function(kappa2, spde, phi, Sigma, mu) {
  Qp <- Q_prime(kappa2, spde)
  log_det_Q <- sum(log(diag(chol(Qp,pivot = T))))
  trace_QEww <- sum(colSums(Qp*Sigma)) + crossprod(mu,Qp)%*%mu
  out <- (trace_QEww / (4*pi*phi) - log_det_Q)@x
  return(out)
}
