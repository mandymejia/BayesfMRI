#' Internal function used in joint approach to group-analysis for combining across models
#'
#' @param theta A vector of hyperparameter values at which to compute the posterior log density
#' @param spde A SPDE object from inla.spde2.matern() function, determines prior precision matrix
#' @param mu_theta Posterior mean from combined subject-level models.
#' @param Q_theta Posterior precision matrix from combined subject-level models.
#' @param M Number of subjects
#' @return A list containing...
#'
#' @inheritSection INLA_Description INLA Requirement
#'
# F.logwt <- function(theta, spde, mu_theta, Q_theta, M){
#   #mu_theta - posterior mean from combined subject-level models
#   #Q_theta - posterior precision matrix from combined subject-level models
#   #M - number of subjects
#   a <- 1; b <- 5e-5
#   n.spde <- (length(theta) - 1)/2
#   mu.tmp <- spde$f$hyper$theta1$param[1:2]
#   mu <- rep(mu.tmp, n.spde)
#   Q.tmp <- matrix(spde$f$hyper$theta1$param[-(1:2)], 2, 2, byrow = TRUE)
#   Q <- kronecker(diag(1, n.spde, n.spde), Q.tmp)
#
#   ## Prior density
#   pr.delta <- dgamma(exp(theta[1]), a, b, log = TRUE) #log prior density on residual precision
#   pr.tk <- as.vector(-t(theta[-1] - mu)%*%Q%*%(theta[-1] - mu))/2 + log(det(Q))/2 - dim(Q)[1]*log(2*pi)/2 #joint log prior density on 2K spde parameters
#   pr.theta <- pr.delta + pr.tk
#
#   (1-M)*pr.theta
# }
