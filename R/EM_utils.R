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
#' @importFrom Matrix diag chol bdiag crossprod colSums
#'
#' @return a scalar
#' @keywords internal
neg_kappa_fn <- function(kappa2, spde, phi, Sigma, mu) {
  Qt <- Q_prime(kappa2, spde)
  KK <- nrow(Sigma) / spde$mesh$n
  if(nrow(Sigma) != spde$mesh$n & KK %% 1 == 0) Qt <- Matrix::bdiag(rep(list(Qt),KK))
  log_det_Q <- sum(2*log(Matrix::diag(Matrix::chol(Qt,pivot = T)))) # compare direct determinant here
  trace_QEww <- (sum(Matrix::colSums(Qt*Sigma)) + Matrix::crossprod(mu,Qt%*%mu))@x
  out <- trace_QEww / (4*pi*phi) - log_det_Q
  return(out)
}

#' Expected log-likelihood function
#'
#' @param Q precision matrix
#' @param sigma2 noise variance
#' @param model_data list with X and y
#' @param Psi data locations to mesh locations matrix
#' @param mu posterior mean of w
#' @param Sigma posterior covariance of w
#' @param A crossprod(X%*%Psi)
#'
#' @return scalar expected log-likelihood
#' @keywords internal
ELL <- function(Q, sigma2, model_data, Psi, mu, Sigma, A) {
  TN <- length(model_data$y)
  R1 <- -TN * log(sigma2) / 2 - crossprod(model_data$y)/(2*sigma2) +
    crossprod(model_data$y,model_data$X%*%Psi%*%mu) -
    Matrix::crossprod(model_data$X%*%Psi%*%mu) / (2*sigma2) -
    sum(Matrix::colSums(A*Sigma)) / (2*sigma2)
  R2 <- sum(2*log(Matrix::diag(Matrix::chol(Q,pivot = T))))/2 -
    sum(Matrix::colSums(Q*Sigma)) / 2 - crossprod(mu,Q%*%mu) / 2
  ELL_out <- R1@x + R2@x
  return(ELL_out)
}
