#' Corrected AIC
#' 
#' Computes corrected AIC (AICc).
#' 
#' @param y The autocorrelated data
#' @param demean Demean \code{y}? Default: \code{FALSE}.
#' @param order.max The model order limit. Default: \code{10}.
#' 
#' @return The cAIC
#' @keywords internal
#' 
#' @importFrom fMRItools is_posNum
AICc <- function(y, demean=FALSE, order.max = 10) {
  N <- length(y)
  stopifnot(is_posNum(order.max))

  # Get regular AIC values.
  ar_mdl <- ar.yw(x = y, aic = FALSE, demean = demean, order.max = order.max)
  AIC_vals <- ar_mdl$aic

  # Get corrected AIC values.
  kseq <- seq(0, order.max)
  AICc <- AIC_vals - (2*(kseq+1)) + 2*N*(kseq+1)/(N-kseq-2)
  AICc <- AICc - min(AICc)

  # Format and return.
  names(AICc) <- paste0("AR(", kseq, ")")
  AICc
}