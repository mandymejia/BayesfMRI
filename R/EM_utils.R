#' Fixed point function for the joint BayesGLMEM update algorithm
#'
#' @param theta a vector containing kappa2, phi, and sigma2, in that order
#' @param spde the spde object
#' @param model_data the model_data object containing \code{y} and \code{X}
#' @param Psi a conversion matrix (N by V) (or N by n)
#' @param K number of covariates
#' @param A The value for Matrix::crossprod(X%*%Psi) (saves time on computation)
#'
#' @return a vector with the same length as \code{theta}, the EM updates
#' @keywords internal
GLMEM_fixptjoint <- function(theta, spde, model_data, Psi, K, A) {
  Q_k <- mapply(BayesfMRI:::spde_Q_phi, kappa2 = theta[1], phi = theta[2], MoreArgs = list(spde = spde), SIMPLIFY = F)
  Q_new <- Matrix::bdiag(rep(Q_k,K))
  Sig_inv <- Q_new + A/theta[3]
  m <- Matrix::crossprod(model_data$X%*%Psi,model_data$y) / theta[3]
  mu <- INLA::inla.qsolve(Sig_inv,m, method = "solve")
  Sigma_new <- INLA::inla.qinv(Sig_inv) # This seems right
  # >> Update sigma2 ----
  sigma2_new <-
    as.numeric(crossprod(model_data$y) -
                 2*Matrix::crossprod(model_data$y,model_data$X%*%Psi%*%mu) +
                 Matrix::crossprod(model_data$X%*%Psi%*%mu) +
                 sum(Matrix::colSums(A*Sigma_new))) / length(model_data$y)
  # >> Update kappa2 ----
  optim_output_k <-
    optimize(
      f = BayesfMRI:::neg_kappa_fn,
      spde = spde,
      phi = theta[2],
      Sigma = Sigma_new,
      mu = mu,
      lower = 0,
      upper = 50
    )
  kappa2_new <- optim_output_k$minimum
  # >> Update phi ----
  Q_tildek <- BayesfMRI:::Q_prime(kappa2_new, spde)
  Q_tilde <- Matrix::bdiag(rep(list(Q_tildek),K))
  Tr_QEww <- (sum(Matrix::colSums(Q_tilde*Sigma_new)) +
                Matrix::crossprod(mu,Q_tilde%*%mu))@x
  phi_new <- Tr_QEww / (4*pi*spde$n.spde*K)
  return(c(kappa2_new,phi_new,sigma2_new))
}

EM_update <-
  function(theta,
           spde,
           Phi,
           K,
           A,
           model_data,
           pct_tol,
           EM_method,
           use_SQUAREM) {
    if(EM_method == "joint") em_fn <- GLMEM_fixptjoint
    if(EM_method == "separate") em_fn <- GLMEM_fixptseparate
    if(use_SQUAREM) {
      squareem_output <-
        squarem(
          par = theta,
          fixptfn = GLMEM_fixptjoint,
          objfn = GLMEM_objfn,
          control = list(tol = tol, trace = verbose),
          spde = spde,
          model_data = model_data,
          Psi = Psi,
          K = K,
          A = A
        )
      theta_new <- squareem_output$par
      kappa2_new <- theta_new[1]
      phi_new <- theta_new[2]
      sigma2_new <- theta_new[3]
    } else {
      step <- 1
      max_pct_change <- Inf
      while(max_pct_change > pct_change_limit) {
        theta_new <-
          GLMEM_fixptjoint(
            theta = theta,
            spde = spde,
            model_data = model_data,
            Psi = Psi,
            K = K,
            A = A
          )
        kappa2_new <- theta_new[1]
        phi_new <- theta_new[2]
        sigma2_new <- theta_new[3]
        sigma2_pct_change <- 100*abs((sigma2_new - sigma2) / sigma2)
        phi_pct_change <- 100*abs((phi_new - phi) / phi)
        kappa2_pct_change <- 100*abs((kappa2_new - kappa2) / kappa2)
        max_pct_change <- max(sigma2_pct_change,phi_pct_change,kappa2_pct_change)
        if(verbose) {
          cat("Step",step, "kappa^2 (%change) =",kappa2_new,"(",kappa2_pct_change,")", "phi (%change) =", phi_new, "(", phi_pct_change,")", "sigma^2 (%change) =",sigma2_new,"(",sigma2_pct_change,")","\n")
        }
        kappa2 <- kappa2_new
        phi <- phi_new
        sigma2 <- sigma2_new
        theta <- theta_new
        step <- step+1
      }
    }
}

#' Fixed point function for the joint BayesGLMEM update algorithm
#'
#' @param theta a list containing kappa2, phi, and sigma2, in that order
#' @param spde the spde object
#' @param model_data the model_data object containing \code{y} and \code{X}
#' @param Psi a conversion matrix (N by V) (or N by n)
#' @param K number of covariates
#' @param A The value for Matrix::crossprod(X%*%Psi) (saves time on computation)
#' @param num.threads (optional) allows users to specify the number of threads used
#'   to work in parallel across the different tasks
#'
#' @importFrom Matrix bdiag colSums crossprod
#' @importFrom INLA inla.qinv inla.qsolve
#' @importFrom parallel detectCores makeCluster parSapply
#'
#' @return a vector with the same length as \code{theta}, the EM updates
#' @keywords internal
GLMEM_fixptseparate <- function(theta, spde, model_data, Psi, K, A, num.threads = 1) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The separate update requires the `parallel` package. Please install it.", call. = FALSE)
  }
  max_num.threads <- min(parallel::detectCores() - 1, 25)
  num.threads <- min(max_num.threads, num.threads)
  cl <- parallel::makeCluster(num.threads)
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
  Sig_inv <- Q_new + A/theta[sigma2_ind]
  m <- Matrix::crossprod(model_data$X%*%Psi,model_data$y) / theta[sigma2_ind]
  mu <- INLA::inla.qsolve(Q = Sig_inv,B = m, method = "solve")
  Sigma_new <- INLA::inla.qinv(Sig_inv)
  # >> Update sigma2 ----
  sigma2_new <-
    as.numeric(crossprod(model_data$y) -
                 2*Matrix::crossprod(model_data$y,model_data$X%*%Psi%*%mu) +
                 Matrix::crossprod(model_data$X%*%Psi%*%mu) +
                 sum(Matrix::colSums(A*Sigma_new))) / length(model_data$y)
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
                   Sigma_new,
                   mu) {
      k_inds <- seq(spde$n.spde) + (k - 1) * spde$n.spde
      # >> Update kappa2 ----
      optim_output_k <-
        optimize(
          f = neg_kappa_fn,
          spde = spde,
          phi = theta[phi_inds][k],
          Sigma = Sigma_new[k_inds, k_inds],
          mu = mu[k_inds, ],
          lower = 0,
          upper = 50
        )
      kappa2_new <- optim_output_k$minimum
      # >> Update phi ----
      Q_tildek <-
        Q_prime(kappa2_new, spde)
      Tr_QEww <-
        (sum(Matrix::colSums(Q_tildek * Sigma_new[k_inds, k_inds])) +
           Matrix::crossprod(mu[k_inds, ], Q_tildek %*%
                               mu[k_inds, ]))@x
      phi_new <-
        Tr_QEww / (4 * pi * spde$n.spde)
      return(c(kappa2_new, phi_new))
    },
    spde = spde,
    theta = theta,
    kappa2_inds = kappa2_inds,
    phi_inds = phi_inds,
    Sigma_new = Sigma_new,
    mu = mu
  )
  parallel::stopCluster(cl)
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
#'
#' @return A scalar value for the negative expected log-likelihood
#' @keywords internal
GLMEM_objfn <- function(theta, spde, model_data, Psi, K, A, num.threads = NULL) {
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
  Q_k <- mapply(BayesfMRI:::spde_Q_phi,
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
  mu <- INLA::inla.qsolve(Sig_inv,m, method = "solve")
  Sigma <- INLA::inla.qinv(Sig_inv)
  R1 <- -TN * log(theta[sigma2_ind]) / 2 - crossprod(model_data$y)/(2*theta[sigma2_ind]) +
    Matrix::crossprod(model_data$y,model_data$X%*%Psi%*%mu) -
    Matrix::crossprod(model_data$X%*%Psi%*%mu) / (2*theta[sigma2_ind]) -
    sum(Matrix::colSums(A*Sigma)) / (2*theta[sigma2_ind])
  R2 <- sum(2*log(Matrix::diag(Matrix::chol(Q,pivot = T))))/2 -
    sum(Matrix::colSums(Q*Sigma)) / 2 - Matrix::crossprod(mu,Q%*%mu) / 2
  ELL_out <- R1@x + R2@x
  return(-ELL_out)
}

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
