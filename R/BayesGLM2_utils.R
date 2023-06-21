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
#' @param nsamp_beta The number of samples to draw from full conditional of beta given the current value of theta (p(beta|theta,y))
#'
#' @importFrom excursions excursions.mc
#' @importFrom Matrix Diagonal
#'
#' @return A list containing \code{mu}, \code{quantiles}, and \code{F}
#'
#' @keywords internal
beta.posterior.thetasamp <- function(
  theta, spde, Xcros, Xycros, contrasts,
  quantiles, excursion_type, gamma, alpha, nsamp_beta=100){

  use_INLA <- TRUE

  n.mesh <- spde$n.spde

  prec.error <- exp(theta[1])
  theta_spde <- matrix(theta[-1], nrow=2) #2xK matrix of the hyperparameters (2 per task)
  K <- ncol(theta_spde)
  M <- length(Xcros)

  use_EM <- all(sapply(c("M0","M1","M2"), function(x) x %in% names(spde)))

  #contruct prior precision matrix for beta, Q_theta for given sampled value of thetas
  # For EM
  if (use_EM) {
    Q.beta <- apply(theta_spde, 2, function(theta_k) {
      theta_k <- exp(theta_k) ^ 2
      out <-
        theta_k[1] * (theta_k[2] ^ 2 * spde$M0 + theta_k[2] * spde$M1 + spde$M2)
      return(out)
    })
  }
  # For INLA
  if (!use_EM) {
    Q.beta <- list()
    for (k in 1:K) {
      theta_k <-
        theta_spde[, k] #theta[(2:3) + 2*(k-1)] #1:2, 2:3, 4:5, ...
      Q.beta[[k]] <-
        INLA::inla.spde2.precision(spde, theta = theta_k) # prior precision for a single task k
    }
  }
  Q <- Matrix::bdiag(Q.beta) #Q_theta in the paper

  N <- dim(Q.beta[[1]])[1] #number of mesh locations
  if(N != n.mesh) stop('Length of betas does not match number of vertices in mesh. Inform developer.')

  beta.samples <- NULL
  # ~5 seconds per subject with PARDISO
  nS <- 1
  Q_nn <- Q
  for(nn in 1:M){
    if(nrow(Q) != nrow(Xcros[[nn]])) {
      nS <- nrow(Xcros[[nn]]) / nrow(Q)
      Q_nn <- Matrix::bdiag(rep(list(Q),nS))
    }
    # compute posterior mean and precision of beta|theta
    Q_nn <- prec.error*Xcros[[nn]] + Q_nn #Q_m in paper
    cholQ_nn <- Matrix::Cholesky(Q_nn)
    if (use_INLA) {
      mu_nn <- INLA::inla.qsolve(Q_nn, prec.error*Xycros[[nn]]) #mu_m in paper
      # draw samples from pi(beta_m|theta,y)
      beta_samp_nn <- INLA::inla.qsample(n = nsamp_beta, Q = Q_nn, mu = mu_nn)
    } else {
      mu_nn <- Matrix::solve(cholQ_nn, prec.error*Xycros[[nn]], system = "A")
      beta_samp_nn <- cholQsample(n = nsamp_beta, cholQ = Q_nn, mu = mu_nn)
    }

    # concatenate samples over models
    beta.samples <- rbind(beta.samples, beta_samp_nn)
  }

  if (excursion_type[1] == 'none') do_excur <- FALSE else do_excur <- TRUE

  # Loop over contrasts
  nC <- length(contrasts)
  mu.contr <- matrix(NA, nrow=n.mesh, ncol=nC)
  if(do_excur) F.contr <- mu.contr else F.contr <- NULL
  if(!is.null(quantiles)){
    num_quantiles <- length(quantiles)
    quantiles.contr <- rep(list(mu.contr), num_quantiles)
    names(quantiles.contr) <- quantiles
  } else {
    num_quantiles <- 0
    quantiles.contr <- NULL
  }

  for (cc in 1:nC) {
    #Construct "A" matrix from paper (linear combinations)
    ctr.mat <- kronecker(t(contrasts[[cc]]), Diagonal(n.mesh, 1))

    #beta.mean.pop.contr <- as.vector(ctr.mat%*%beta.mean.pop.mat)  # NKx1 or Nx1
    samples_cc <- as.matrix(ctr.mat%*%beta.samples)  # N x nsamp_beta
    mu.contr[,cc] <- rowMeans(samples_cc) #compute mean over beta samples
    if(num_quantiles > 0){
      for(iq in 1:num_quantiles){
        quantiles.contr[[iq]][,cc] <- apply(samples_cc, 1, quantile, quantiles[iq])
      }
    }

    # Estimate excursions set for current contrast
    if (do_excur) {
      excur_cc <- excursions::excursions.mc(
        samples_cc, u = gamma[cc], type = excursion_type[cc], alpha = alpha[cc]
      )
      F.contr[,cc] <- excur_cc$F
    }
  }

  list(
    mu = mu.contr,
    quantiles = quantiles.contr,
    F = F.contr
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
  mu.tmp <- spde$f$hyper$theta1$param[1:2]
  mu <- rep(mu.tmp, n.spde)
  Q.tmp <- matrix(spde$f$hyper$theta1$param[-(1:2)], 2, 2, byrow = TRUE)
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

#' Summarize a \code{"BayesGLM2"} object
#'
#' Summary method for class \code{"BayesGLM2"}
#'
#' @param object Object of class \code{"BayesGLM2"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BayesGLM2"} object, a list summarizing the 
#'  properties of \code{object}.
#' @method summary BayesGLM2
summary.BayesGLM2 <- function(object, ...) {

  x <- list(
    n_contrasts = length(object$contrasts),
    tasks = object$task_names,
    sessions = object$session_names,
    n_loc_total = vapply(lapply(object$model_results, '[[', "mask"), length, 0),
    n_loc_modeled = vapply(lapply(object$model_results, '[[', "mask"), sum, 0),
    excursion_type = object$excursion_type
  )
  class(x) <- "summary.BayesGLM2"
  return(x)
}

#' @rdname summary.BayesGLM2
#' @export
#'
#' @param x Object of class \code{"summary.BayesGLM2"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BayesGLM2
print.summary.BayesGLM2 <- function(x, ...) {
  cat("====BayesGLM2 result===================\n")
  cat("Tasks:    ", paste0("(", length(x$tasks), ") ", paste(x$tasks, collapse=", ")), "\n")
  if (length(x$sessions)==1 && x$sessions == "session_combined") {
    cat("Sessions: ", paste0("(", x$n_sess_orig, ", combined) \n"))
  } else {
    cat("Sessions: ", paste0("(", length(x$sessions), ") ", paste(x$sessions, collapse=", ")), "\n")
  }
  cat("Locations:\n")
  for (ii in seq(length(x$n_loc_total))) {
    cat(
      "          ", paste0(names(x$n_loc_total)[ii], ": ", x$n_loc_modeled[[ii]]),
      "modeled,", x$n_loc_total[[ii]], "total", "\n"
    )
  }
  cat("Excursion:", paste(x$excursion_type, collapse=", "), "\n")
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.BayesGLM2
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BayesGLM2
print.BayesGLM2 <- function(x, ...) {
  print.summary.BayesGLM2(summary(x))
}

#' Summarize a \code{"BayesGLM2_cifti"} object
#'
#' Summary method for class \code{"BayesGLM2_cifti"}
#'
#' @param object Object of class \code{"BayesGLM2_cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BayesGLM2_cifti"} object, a list summarizing the 
#'  properties of \code{object}.
#' @method summary BayesGLM2_cifti
summary.BayesGLM2_cifti <- function(object, ...) {
  x <- summary.BayesGLM2(object$BayesGLM2_results)
  class(x) <- "summary.BayesGLM2_cifti"
  x
}

#' @rdname summary.BayesGLM2_cifti
#' @export
#'
#' @param x Object of class \code{"summary.BayesGLM2_cifti"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BayesGLM2_cifti
print.summary.BayesGLM2_cifti <- function(x, ...) {
  print.summary.BayesGLM2(x, ...)
}

#' @rdname summary.BayesGLM2_cifti
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BayesGLM2_cifti
print.BayesGLM2_cifti <- function(x, ...) {
  print.summary.BayesGLM2_cifti(summary(x))
}
