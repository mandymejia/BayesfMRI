
#' Summarise BayesGLM objects
#'
#' Summary method for class "BayesGLM"
#'
#' @param object an object of class "BayesGLM"
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary BayesGLM
summary.BayesGLM <- function(object, ...)
{
  out <- list()
  class(out) <- "summary.BayesGLM"
  out$sessions <- object$session_names
  out$betas <- object$beta_names
  out$call <- object$INLA_result$call
  out$inla.summary <- summary(object$model)
  return(out)
}


#' @param x an object of class "summary.BayesGLM"
#' @export
#' @method print summary.BayesGLM
#' @rdname summary.BayesGLM
print.summary.BayesGLM <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("Sessions: ", x$sessions,"\n")
  cat("Time used:\n", x$inla.summary$cpu.used)
}

#' @export
#' @method print BayesGLM
#' @rdname summary.BayesGLM
print.BayesGLM <- function(x, ...) {
  print.summary.BayesGLM(summary(x))
}

# TO DO: Add print and summary functions for session object (may need as.session function too, and move is.session here)



#' Find nonzero element in a matrix using 2-means clustering
#'
#' @param beta_est A vector or matrix of values from which values close to zero should be assigned a value of zero.
#'
#' @return A vector or matrix of the same dimension as beta_est in which values close to zero are assigned the value of zero. The closeness of a value to zero is found by performing two-means clustering on the absolute values of beta_est, and
#' @export
#'
find_nonzero <- function(beta_est) {
  vector_beta <- c(beta_est)
  if(any(is.na(vector_beta))) vector_beta <- vector_beta[!is.na(vector_beta)]
  km_beta <- kmeans(abs(vector_beta),2)
  which_nonzero <- which.max(km_beta$centers[,1])
  keep_nonzero <- as.numeric(km_beta$cluster == which_nonzero)
  out <- beta_est
  out[!is.na(out)] <- out[!is.na(out)] * keep_nonzero
  return(out)
}

#' Sequential 2-means variable selection
#'
#' @param x A vector consisting of all variables of interest for a single draw from a posterior distribution
#' @param b A scale parameter used to determine at what distance cluster centers are considered to be the same.
#'
#' @return The number of nonzero values detected within x
#' @export
#'
s2m <- function(x,b){
  two_means <- kmeans(abs(x),2)
  zero_idx <- which(two_means$cluster == which.min(two_means$centers))
  A <- x[zero_idx]
  two_centers <- kmeans(abs(A),2,algorithm=c("Lloyd"))
  iterations <- 1
  while(abs(two_centers$centers[1, 1] - two_centers$centers[2, 1]) > b) {
    zero_idx <- which(two_centers$cluster == which.min(two_centers$centers))
    A <- A[zero_idx]
    two_centers <- kmeans(abs(A),2,algorithm=c("Lloyd"))
    iterations <- iterations + 1
  }
  num_nonzero <- length(x) - length(A)
  return(num_nonzero)
}

#' Sequential 2-means on array B
#'
#' @param B An array of posterior samples (typically a matrix), in which the last margin corresponds to a single posterior sample
#' @param sigma  A scale parameter used to determine at what distance cluster centers are considered to be the same.
#'
#' @return An array of dimension `head(dim(B),-1)` with a point estimate of B based on the sequential 2-means method
#' @export
#'
#' @md
s2m_B <- function(B,sigma){
  nonzero_nums <- sapply(asplit(B,length(dim(B))),function(B_s) s2m(c(B_s),sigma))
  num_nonzero <- ceiling(median(nonzero_nums))
  median_B <- apply(B,seq(length(dim(B)) - 1),median)
  cutoff <- quantile(c(abs(median_B)),1 - num_nonzero/length(median_B))
  out <- median_B
  out[which(out < cutoff)] <- 0
  return(out)
}


