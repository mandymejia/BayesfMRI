#' Beta posterior theta sampling
#'
#' Internal function used in joint approach to group-analysis
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param theta A single sample of theta (hyperparameters) from q(theta|y)
#' @param spde A SPDE object from inla.spde2.matern() function.
#' @param Xcros A crossproduct of design matrix.
#' @param Xycros A crossproduct of design matrix and BOLD y.
#' @param contrasts A list of vectors of length M*K specifying the contrasts of interest.
#' @param quantiles Vector of posterior quantiles to return in addition to the posterior mean
#' @param excursion_type Vector of excursion function type (">", "<", "!=") for each contrast
#' @param gamma Vector of activation thresholds for each contrast
#' @param alpha Significance level for activation for the excursion sets
#' @param alpha_grid Optional list of numeric vectors specifying multiple alpha thresholds for each contrast for nested excursion set inference. If provided, `alpha` is ignored.
#' @param nsamp_beta The number of samples to draw from full conditional of beta given the current value of theta (p(beta|theta,y))
#'
#' @importFrom excursions excursions.mc
#' @importFrom Matrix Diagonal
#'
#' @return A list containing \code{mu}, \code{quantiles}, and \code{F}
#'
#' @keywords internal
beta.posterior.thetasamp <- function(
    theta,
    spde,
    Xcros,
    Xycros,
    contrasts,
    quantiles,
    excursion_type,
    gamma,
    alpha = NULL,
    alpha_grid = NULL,
    nsamp_beta = 100
){

  n.mesh <- spde$n.spde

  prec.error <- exp(theta[1])
  theta_spde <- matrix(theta[-1], nrow = 2)  # 2 x K matrix of hyperparameters
  K <- ncol(theta_spde)
  M <- length(Xcros)

  use_EM <- all(sapply(c("M0", "M1", "M2"), function(x) x %in% names(spde)))

  # Construct prior precision matrix for beta, Q_theta for given sampled value of theta
  # For EM
  if (use_EM) {
    Q.beta <- apply(theta_spde, 2, function(theta_k) {
      theta_k <- exp(theta_k)^2
      out <- theta_k[1] * (theta_k[2]^2 * spde$M0 + theta_k[2] * spde$M1 + spde$M2)
      return(out)
    })
  }

  # For INLA
  if (!use_EM) {
    Q.beta <- list()
    for (k in 1:K) {
      theta_k <- theta_spde[, k]
      Q.beta[[k]] <- INLA::inla.spde2.precision(spde, theta = theta_k)
    }
  }

  N <- dim(Q.beta[[1]])[1]  # number of mesh locations
  if (N != n.mesh) {
    stop("Length of betas does not match number of vertices in mesh. Inform developer.")
  }

  beta.samples <- NULL

  # ~5 seconds per subject with PARDISO
  nS <- 1
  Q <- Q_theta <- Matrix::bdiag(Q.beta)  # Q_theta in the paper

  for (mm in seq(M)) {
    if (nrow(Q) != nrow(Xcros[[mm]])) {
      nS <- nrow(Xcros[[mm]]) / nrow(Q)
      if (nS != round(nS)) {
        stop("Internal error.")
      }
      Q_theta <- Matrix::bdiag(rep(list(Q), nS))
    }

    # Compute posterior mean and precision of beta | theta
    Q_mm <- prec.error * Xcros[[mm]] + Q_theta
    cholQ_mm <- Matrix::Cholesky(Q_mm)
    mu_mm <- INLA::inla.qsolve(Q_mm, prec.error * Xycros[[mm]])

    # Draw samples from pi(beta_m | theta, y)
    beta_samp_mm <- INLA::inla.qsample(n = nsamp_beta, Q = Q_mm, mu = mu_mm)

    # Concatenate samples over models
    beta.samples <- rbind(beta.samples, beta_samp_mm)
  }

  do_excur <- !(excursion_type[1] == "none")

  # Loop over contrasts
  nC <- length(contrasts)
  mu.contr <- matrix(NA, nrow = n.mesh, ncol = nC)

  use_nested_alpha <- !is.null(alpha_grid)

  # Validate alpha / alpha_grid only when excursion inference is requested
  if (do_excur) {
    if (length(gamma) == 1) {
      gamma <- rep(gamma, nC)
    }
    if (length(gamma) != nC) {
      stop("Length of `gamma` must match number of contrasts or be equal to one.")
    }

    if (length(excursion_type) == 1) {
      excursion_type <- rep(excursion_type, nC)
    }
    if (length(excursion_type) != nC) {
      stop("Length of `excursion_type` must match number of contrasts or be equal to one.")
    }

    if (use_nested_alpha) {
      # alpha_grid can be:
      #   (1) one numeric vector shared by all contrasts
      #   (2) a list of numeric vectors, one per contrast
      if (is.numeric(alpha_grid)) {
        alpha_grid <- replicate(
          nC,
          sort(unique(alpha_grid), decreasing = TRUE),
          simplify = FALSE
        )
      } else if (is.list(alpha_grid)) {
        if (length(alpha_grid) != nC) {
          stop("`alpha_grid` must have length equal to the number of contrasts.")
        }
        alpha_grid <- lapply(alpha_grid, function(x) {
          if (!is.numeric(x)) {
            stop("Each element of `alpha_grid` must be numeric.")
          }
          sort(unique(x), decreasing = TRUE)
        })
      } else {
        stop("`alpha_grid` must be either a numeric vector or a list of numeric vectors.")
      }

      bad_len <- vapply(alpha_grid, length, integer(1)) == 0
      if (any(bad_len)) {
        stop("Each contrast must have at least one alpha threshold in `alpha_grid`.")
      }

      bad_val <- vapply(
        alpha_grid,
        function(x) any(!is.finite(x) | x <= 0 | x >= 1),
        logical(1)
      )
      if (any(bad_val)) {
        stop("All values in `alpha_grid` must be finite and strictly between 0 and 1.")
      }

      # In nested mode, alpha is ignored
      alpha <- NULL

    } else {
      # Single-threshold mode
      if (is.null(alpha)) {
        stop("`alpha` must be provided in single-threshold mode.")
      }
      if (length(alpha) == 1) {
        alpha <- rep(alpha, nC)
      }
      if (length(alpha) != nC) {
        stop("Length of `alpha` must match number of contrasts or be equal to one.")
      }
      if (any(!is.finite(alpha) | alpha <= 0 | alpha >= 1)) {
        stop("All values in `alpha` must be finite and strictly between 0 and 1.")
      }
    }
  }

  # Initialize outputs
  if (do_excur) {
    if (use_nested_alpha) {
      F.contr <- NULL
      F.contr_levels <- vector("list", nC)
      names(F.contr_levels) <- names(contrasts)
    } else {
      F.contr <- mu.contr
      F.contr_levels <- NULL
    }
  } else {
    F.contr <- NULL
    F.contr_levels <- NULL
  }

  if (!is.null(quantiles)) {
    num_quantiles <- length(quantiles)
    quantiles.contr <- rep(list(mu.contr), num_quantiles)
    names(quantiles.contr) <- quantiles
  } else {
    num_quantiles <- 0
    quantiles.contr <- NULL
  }

  for (cc in seq_len(nC)) {

    # Construct "A" matrix from paper (linear combinations)
    ctr.mat <- kronecker(t(contrasts[[cc]]), Matrix::Diagonal(n.mesh, 1))

    # N x nsamp_beta
    samples_cc <- as.matrix(ctr.mat %*% beta.samples)

    # Posterior mean over beta samples
    mu.contr[, cc] <- rowMeans(samples_cc)

    # Posterior quantiles over beta samples
    if (num_quantiles > 0) {
      for (iq in seq_len(num_quantiles)) {
        quantiles.contr[[iq]][, cc] <- apply(samples_cc, 1, quantile, quantiles[iq])
      }
    }

    # Excursion calculations
    if (do_excur) {
      if (use_nested_alpha) {

        aa <- alpha_grid[[cc]]

        F_cc_levels <- vapply(
          aa,
          function(a) {
            excur_tmp <- excursions::excursions.mc(
              samples_cc,
              u = gamma[cc],
              type = excursion_type[cc],
              alpha = a
            )
            excur_tmp$F
          },
          numeric(n.mesh)
        )

        colnames(F_cc_levels) <- paste0("alpha_", aa)
        F.contr_levels[[cc]] <- F_cc_levels

      } else {

        excur_cc <- excursions::excursions.mc(
          samples_cc,
          u = gamma[cc],
          type = excursion_type[cc],
          alpha = alpha[cc]
        )
        F.contr[, cc] <- excur_cc$F
      }
    }
  }

  list(
    mu = mu.contr,
    quantiles = quantiles.contr,
    F = F.contr,
    F_levels = F.contr_levels
  )
}


#' F logwt
#'
#' Internal function used in joint approach to group-analysis for combining across models
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param theta A vector of hyperparameter values at which to compute the posterior log density
#' @param spde A SPDE object from inla.spde2.matern() function, determines prior precision matrix
#' @param mu_theta Posterior mean from combined subject-level models.
#' @param Q_theta Posterior precision matrix from combined subject-level models.
#' @param nN Number of subjects
#' @return The prior density
#'
#' @importFrom stats dgamma
#'
#' @keywords internal
F.logwt <- function(theta, spde, mu_theta, Q_theta, nN){
  #mu_theta - posterior mean from combined subject-level models
  #Q_theta - posterior precision matrix from combined subject-level models
  #nN - number of subjects
  a <- 1; b <- 5e-5
  n.spde <- (length(theta) - 1)/2
  mu.tmp <- spde$f$hyper$theta1$param[1:2] #prior mean for (log(tau), log(kappa))
  mu <- rep(mu.tmp, n.spde)
  Q.tmp <- matrix(spde$f$hyper$theta1$param[-(1:2)], 2, 2, byrow = TRUE) #prior precision matrix for (log(tau), log(kappa))
  Q <- kronecker(diag(1, n.spde, n.spde), Q.tmp)

  ## Prior density
  pr.delta <- dgamma(exp(theta[1]), a, b, log = TRUE) #log prior density on residual precision
  pr.tk <- as.vector(-t(theta[-1] - mu)%*%Q%*%(theta[-1] - mu))/2 + log(det(Q))/2 - dim(Q)[1]*log(2*pi)/2 #joint log prior density on 2K spde parameters
  pr.theta <- pr.delta + pr.tk

  (1-nN)*pr.theta
}

#' Sample from a multivariate normal with mean and precision
#'
#' @param n number of samples
#' @param mu mean vector (length = p)
#' @param Q sparse p x p positive definite precision matrix (class = dgCMatrix)
#'
#' @return An n x p matrix of samples
#'
#' @importFrom Matrix solve
#' @importFrom stats rnorm
#' @keywords internal
qsample <- function(n, mu, Q) {
  p <- length(mu)
  if(p != nrow(Q) | p != ncol(Q)) stop("Dimension mismatch between mu and Q.")
  cholQ <- Matrix::Cholesky(Q)
  Z <- matrix(rnorm(n*p), nrow = n, ncol = p)
  out <- Matrix::solve(cholQ,Z, system = "A")
  out <- out + mu
  return(out)
}

#' Sample from the multivariate normal distribution with Cholesky(Q)
#'
#' @param n number of samples
#' @param mu mean vector
#' @param cholQ Cholesky decomposition of the precision (found via \code{Matrix::Cholesky(Q)})
#'
#' @return An \eqn{n \times p} matrix of samples from the MVN distribution,
#'  where \eqn{p} is the length of \code{mu}.
#'
#' @importFrom stats rnorm
#' @importFrom Matrix solve
#' @keywords internal
cholQsample <- function(n, mu, cholQ) {
  p <- length(mu)
  if(p != nrow(cholQ) | p != ncol(cholQ)) stop("Dimension mismatch between mu and Q.")
  Z <- matrix(rnorm(n*p), nrow = n, ncol = p)
  out <- Matrix::solve(cholQ,Z, system = "A")
  out <- out + mu
  return(out)
}
