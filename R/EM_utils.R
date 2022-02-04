#' Fixed point function for the joint BayesGLMEM update algorithm
#'
#' @param theta a list containing kappa2, phi, and sigma2, in that order
#' @param spde the spde object
#' @param model_data the model_data object containing \code{y} and \code{X}
#' @param Psi a conversion matrix (N by V) (or N by n)
#' @param K number of covariates
#' @param A The value for Matrix::crossprod(X%*%Psi) (saves time on computation)
#' @param cl parallelization cluster
#' @param Ns The number of samples used to approximate traces using the Hutchinson
#'   estimator. If set to 0, the exact trace is found.
#'
#' @importFrom Matrix bdiag colSums crossprod solve
#' @importFrom parallel detectCores makeCluster parSapply
#' @importFrom stats optimize
#'
#' @return a vector with the same length as \code{theta}, the EM updates
#' @keywords internal
GLMEM_fixptseparate <- function(theta, spde, model_data, Psi, K, A, cl, Ns = 50) {
  kappa2_inds <- seq(K)
  phi_inds <- seq(K) + K
  sigma2_ind <- 2*K + 1
  Q_k <-
    mapply(
      spde_Q_phi,
      kappa2 = theta[kappa2_inds],
      phi = theta[phi_inds],
      MoreArgs = list(spde = spde),
      SIMPLIFY = F
    )
  Q_new <- Matrix::bdiag(Q_k)
  n_sess_em <- nrow(A) / nrow(Q_new)
  if(n_sess_em > 1) Q_new <- Matrix::bdiag(lapply(seq(n_sess_em),function(x) Q_new))
  Sig_inv <- Q_new + A/theta[sigma2_ind]
  m <- Matrix::crossprod(model_data$X%*%Psi,model_data$y) / theta[sigma2_ind]
  mu <- Matrix::solve(Sig_inv, m)
  X_Psi_mu <- model_data$X%*%Psi%*%mu
  cp_X_Psi_mu <- Matrix::crossprod(X_Psi_mu)
  # >> Update sigma2 ----
  if(Ns == 0) {
    Sigma_new <- Matrix::solve(Sig_inv)
    traceAEww <-
      cp_X_Psi_mu +
      sum(Matrix::colSums(A*Sigma_new))
    kappa_fn <- neg_kappa_fn
  }
  if(Ns > 0) {
    Vh <- matrix(sample(x = c(-1,1), size = Ns * nrow(A), replace = TRUE),
                 nrow(A), Ns)
    P <- Matrix::solve(Sig_inv, Vh)
    traceAEww <-
      as.numeric(Matrix::crossprod(mu,A %*% mu)) +
      TrSigB(P,A,Vh)
    kappa_fn <- neg_kappa_fn2
  }
  sigma2_new <-
    as.numeric(crossprod(model_data$y) -
                 2*Matrix::crossprod(model_data$y,X_Psi_mu) +
                 traceAEww
    ) / length(model_data$y)
  kappa2_new <- theta[kappa2_inds]
  phi_new <- theta[phi_inds]
  kp <- parallel::parSapply(
    cl = cl,
    seq(K),
    FUN = function(k,
                   spde,
                   theta,
                   kappa2_inds,
                   phi_inds,
                   P,
                   mu,
                   Vh
    ) {
      # source("~/github/BayesfMRI/R/EM_utils.R") # For debugging
      big_K <- length(kappa2_inds)
      big_N <- spde$n.spde
      n_sess_em <- length(mu) / (big_K * big_N)
      k_inds <- c(sapply(seq(n_sess_em), function(ns) {
        seq( big_N * (big_K * (ns - 1) + k - 1) + 1, big_N * (big_K * (ns - 1) + k))
      }))
      # >> Update kappa2 ----
      optim_output_k <-
        stats::optimize(
          f = neg_kappa_fn2,
          spde = spde,
          phi = theta[phi_inds][k],
          P = P[k_inds,],
          mu = mu[k_inds, ],
          Vh = Vh[k_inds,], # Comment this out if using neg_kappa_fn
          lower = 0,
          upper = 50
        )
      kappa2_new <- optim_output_k$minimum
      # >> Update phi ----
      Tr_QEww <-
        TrQEww(kappa2 = kappa2_new, spde = spde, P = P[k_inds,], mu = mu[k_inds,],Vh = Vh[k_inds,])
      phi_new <-
        Tr_QEww / (4 * pi * spde$n.spde * n_sess_em)
      return(c(kappa2_new, phi_new))
    },
    spde = spde,
    theta = theta,
    kappa2_inds = kappa2_inds,
    phi_inds = phi_inds,
    P = P,
    mu = mu,
    Vh = Vh
  )
  # parallel::stopCluster(cl)
  return(c(kp[1,],kp[2,],sigma2_new))
}

#' Objective function for the BayesGLM EM algorithm
#'
#' This returns the negative of the expected log-likelihood function.
#'
#' @param theta a vector containing kappa2, phi, and sigma2, in that order
#' @param spde the spde object
#' @param model_data the model_data object containing \code{y} and \code{X}
#' @param Psi a conversion matrix (N by V) (or N by n)
#' @param K number of covariates
#' @param A The value for Matrix::crossprod(X%*%Psi) (saves time on computation)
#' @param num.threads Needed for SQUAREM (it is an argument to the fixed-point functions)
#' @param Ns The number of samples used to approximate traces using the Hutchinson
#'   estimator. If set to 0, the exact trace is found.
#'
#' @return A scalar value for the negative expected log-likelihood
#' @keywords internal
GLMEM_objfn <- function(theta, spde, model_data, Psi, K, A, num.threads = NULL, Ns = NULL) {
  if(length(theta) > 3) { # This condition means that parameters are being updated separately
    kappa2_inds <- seq(K)
    phi_inds <- seq(K) + K
    sigma2_ind <- 2 *K + 1
  } else {
    kappa2_inds <- 1
    phi_inds <- 2
    sigma2_ind <- 3
  }
  TN <- length(model_data$y)
  Q_k <- mapply(spde_Q_phi,
                kappa2 = theta[kappa2_inds],
                phi = theta[phi_inds],
                MoreArgs = list(spde = spde),
                SIMPLIFY = F)
  if(length(Q_k) > 1) {
    Q <- Matrix::bdiag(Q_k)
  } else {
    Q <- Matrix::bdiag(rep(Q_k,K))
  }
  Sig_inv <- Q + A/theta[sigma2_ind]
  m <- Matrix::crossprod(model_data$X%*%Psi,model_data$y) / theta[sigma2_ind]
  mu <- Matrix::solve(Sig_inv,m)
  XB <- model_data$X%*%Psi%*%mu
  y_res <- model_data$y - XB
  ELL_out <- as.numeric(-TN * log(theta[sigma2_ind]) / 2 - Matrix::crossprod(y_res))
  return(-ELL_out)
}

#' Calculate the SPDE covariance
#'
#' @param kappa2 A scalar
#' @param phi A scalar
#' @param spde An object containing the information about the
#'   mesh structure for the SPDE prior
#'
#' @return The SPDE prior matrix
#' @keywords internal
spde_Q_phi <- function(kappa2, phi, spde) {
  # Cmat <- spde$M0
  # Gmat <- (spde$M1 + Matrix::t(spde$M1))/2
  # GtCinvG <- spde$M2
  Q <- (kappa2*spde$M0 + (spde$M1 + Matrix::t(spde$M1)) + spde$M2/kappa2) / (4*pi*phi)
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
  Cmat <- spde$M0
  Gmat <- (spde$M1 + Matrix::t(spde$M1))/2
  GtCinvG <- spde$M2
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

#' The negative of the objective function for kappa - Sig_inv
#'
#' @param kappa2 scalar
#' @param spde an spde object
#' @param phi scalar
#' @param P Matrix of dimension nk by N_s found by \code{solve(Sig_inv,Vh)}
#' @param mu dgeMatrix
#' @importFrom Matrix diag chol bdiag
#'
#' @return a scalar
#' @keywords internal
neg_kappa_fn2 <- function(kappa2, spde, phi, P, mu, Vh) {
  Qt <- Q_prime(kappa2, spde)
  n_spde <- spde$mesh$n
  if(is.null(n_spde)) n_spde <- spde$n.spde
  KK <- nrow(P) / n_spde
  if(nrow(P) != n_spde & KK %% 1 == 0) Qt <- Matrix::bdiag(rep(list(Qt),KK))
  log_det_Q <- sum(2*log(Matrix::diag(Matrix::chol(Qt,pivot = T)))) # compare direct determinant here
  trace_QEww <-
    TrQEww(
      kappa2 = kappa2,
      spde = spde,
      P = P,
      mu = mu,
      Vh = Vh
    )
  out <- trace_QEww / (4*pi*phi) - log_det_Q
  return(out)
}

#' Trace approximation function
#'
#' @param kappa2 a scalar
#' @param spde spde object
#' @param P Matrix of dimension nk by N_s found by \code{solve(Sig_inv,Vh)}
#' @param mu posterior mean
#' @param Vh matrix of random variables with \code{nrow(Sig_inv)} rows and Ns
#'   columns
#' @importFrom Matrix t crossprod bdiag
#'
#' @return a scalar
#' @keywords internal
TrQEww <- function(kappa2, spde, P, mu, Vh){
  Cmat <- spde$M0
  Gmat <- (spde$M1 + Matrix::t(spde$M1))/2
  GtCinvG <- spde$M2
  n_spde <- spde$mesh$n
  if(is.null(n_spde)) n_spde <- spde$n.spde
  n_sess_em <- nrow(P) / n_spde
  if(n_sess_em > 1) {
    Cmat <- Matrix::bdiag(rep(list(Cmat),n_sess_em))
    Gmat <- Matrix::bdiag(rep(list(Gmat),n_sess_em))
    GtCinvG <- Matrix::bdiag(rep(list(GtCinvG),n_sess_em))
    # Cmat <- spam::bdiag(rep(list(Cmat),n_sess_em))
    # Gmat <- spam::bdiag(rep(list(Gmat),n_sess_em))
    # GtCinvG <- spam::bdiag(rep(list(GtCinvG),n_sess_em))
  }
  k2uCu <- kappa2*Matrix::crossprod(mu,Cmat%*%mu)
  two_uGu <- 2*Matrix::crossprod(mu,Gmat)%*%mu # This does not depend on kappa2, but needed for phi
  kneg2uGCGu <- (1/kappa2)*Matrix::crossprod(mu,GtCinvG%*%mu)
  # k2uCu <- kappa2*spam::crossprod(mu,Cmat%*%mu)
  # two_uGu <- 2*spam::crossprod(mu,Gmat)%*%mu # This does not depend on kappa2, but needed for phi
  # kneg2uGCGu <- (1/kappa2)*spam::crossprod(mu,GtCinvG%*%mu)
  k2TrCSig <- kappa2 * TrSigB(P,Cmat,Vh)
  twoTrGSig <- 2 * TrSigB(P,Gmat,Vh)
  kneg2GCGSig <- TrSigB(P,GtCinvG,Vh) / kappa2
  trace_QEww <- as.numeric(k2uCu + kneg2uGCGu + k2TrCSig + twoTrGSig + kneg2GCGSig)
  return(trace_QEww)
}

#' Hutchinson estimator of the trace
#'
#' @param P Matrix of dimension nk by N_s found by \code{solve(Sig_inv,Vh)}
#' @param B Matrix of dimension nk by nk inside the desired trace product
#' @param vh Matrix of dimension nk by N_s in which elements are -1 or 1 with
#'   equal probability.
#'
#' @return a scalar estimate of the trace of \code{Sigma %*% B}
#'
#' @importFrom Matrix diag crossprod
#' @keywords internal
TrSigB <- function(P,B,vh) {
  Ns <- ncol(P)
  out <- sum(Matrix::diag(Matrix::crossprod(P,B %*% vh))) / Ns
  # out <- sum(spam::diag(spam::crossprod(P,B %*% vh))) / Ns
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

#' Trace of Q beta' beta
#'
#' @param kappa2 scalar
#' @param beta_hat a vector
#' @param spde an spde object
#'
#' @return a scalar
#' @keywords internal
TrQbb <- function(kappa2, beta_hat, spde) {
  Qt <- Q_prime(kappa2, spde)
  KK <- length(beta_hat) / spde$mesh$n
  if(length(beta_hat) != spde$mesh$n & KK %% 1 == 0) Qt <- Matrix::bdiag(rep(list(Qt),KK))
  out <- sum(diag(Qt %*% tcrossprod(beta_hat)))
  # The below is an approximation that may be faster in higher dimensions
  # Ns <- 50
  # v <- matrix(sample(x = c(-1,1), size = nrow(Qt) * Ns, replace = TRUE), nrow(Qt), Ns)
  # vQ <- apply(v,2, crossprod, y = Qt)
  # vbb <- apply(v,2,crossprod, y = tcrossprod(beta_hat))
  # out2 <- sum(unlist(mapply(function(q,b) (q %*% b)@x, q = vQ, b = split(vbb, col(vbb))))) / Ns
  return(out)
}

#' Function to optimize over kappa2
#'
#' @param kappa2 scalar
#' @param phi scalar
#' @param spde an spde object
#' @param beta_hat vector
#'
#' @return a scalar
#' @keywords internal
kappa_init_fn <- function(kappa2, phi, spde, beta_hat) {
  Qt <- Q_prime(kappa2, spde)
  n_spde <- spde$mesh$n
  if(is.null(n_spde)) n_spde <- spde$n.spde
  KK <- length(beta_hat) / n_spde
  if(length(beta_hat) != n_spde & KK %% 1 == 0) Qt <- Matrix::bdiag(rep(list(Qt),KK))
  log_det_Q <- sum(2*log(Matrix::diag(Matrix::chol(Qt,pivot = T))))
  # out <- log_det_Q / 2 - TrQbb(kappa2,beta_hat,spde) / (8*pi*phi)
  out <- log_det_Q / 2 - as.numeric(Matrix::crossprod(beta_hat, Qt %*% beta_hat)) / (8*pi*phi)
  return(-out)
}

#' Objective function for the initialization of kappa2 and phi
#'
#' @param theta a vector c(kappa2,phi)
#' @param spde an spde object
#' @param beta_hat vector
#'
#' @return scalar
#' @keywords internal
init_objfn <- function(theta, spde, beta_hat) {
  QQ <- spde_Q_phi(kappa2 = theta[1],phi = theta[2], spde)
  KK <- length(beta_hat) / spde$mesh$n
  if(length(beta_hat) != spde$mesh$n & KK %% 1 == 0) QQ <- Matrix::bdiag(rep(list(QQ),KK))
  log_det_Q <- sum(2*log(Matrix::diag(Matrix::chol(QQ,pivot = T))))
  out <- (log_det_Q / 2 - crossprod(beta_hat,QQ)%*%beta_hat / 2)@x
  return(-out)
}

#' The fix point function for the initialization of kappa2 and phi
#'
#' @param theta a vector c(kappa2,phi)
#' @param spde an spde object
#' @param beta_hat vector
#'
#' @importFrom stats optimize
#'
#' @return scalar
#' @keywords internal
init_fixpt <- function(theta, spde, beta_hat) {
  kappa2 <- theta[1]
  n_spde <- spde$mesh$n
  if(is.null(n_spde)) n_spde <- spde$n.spde
  num_sessions <- length(beta_hat) / n_spde
  Qp <- Q_prime(kappa2, spde)
  if(num_sessions > 1) Qp <- Matrix::bdiag(rep(list(Qp), num_sessions))
  phi <- as.numeric(Matrix::crossprod(beta_hat, Qp %*% beta_hat)) / (4 * pi * n_spde * num_sessions)
  kappa2 <-
    stats::optimize(
      f = kappa_init_fn,
      phi = phi,
      spde = spde,
      beta_hat = beta_hat,
      lower = 0,
      upper = 50,
      maximum = FALSE
    )$minimum
  return(c(kappa2, phi))
}
