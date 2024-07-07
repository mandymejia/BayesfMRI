#' Estimate residual autocorrelation for prewhitening
#'
#' @param resids Estimated residuals in \eqn{T \times V} numeric matrix
#' @param ar_order,aic Order of the AR model used to prewhiten the data at each location.
#'  If \code{!aic} (default), the order will be exactly \code{ar_order}. If \code{aic},
#'  the order will be between zero and \code{ar_order}, as determined by the AIC.
#' @importFrom stats ar.yw
#'
#' @keywords internal
#'
#' @return Estimated AR coefficients and residual variance at every vertex
pw_estimate <- function(resids, ar_order, aic=FALSE){

  V <- ncol(resids)
  AR_coefs <- matrix(NA, V, ar_order)
  AR_resid_var <- rep(NA, V)
  AR_AIC <- if (aic) {rep(NA, V) } else { NULL }
  for (v in seq(V)) {
    if (is.na(resids[1,v])) { next }

    # # If `AIC`, overwrite the model order with the one selected by `cAIC`.
    # if (aic) { ar_order <- which.min(cAIC(resids, order.max=ar_order)) - 1 }

    ar_v <- ar.yw(resids[,v], aic = aic, order.max = ar_order)
    aic_order <- ar_v$order # same as length(ar_v$ar)
    AR_coefs[v,] <- c(ar_v$ar, rep(0, ar_order-aic_order)) # The AR parameter estimates
    AR_resid_var[v] <- ar_v$var.pred # Residual variance
    if (aic) { AR_AIC[v] <- ar_v$order } # Model order
  }

  list(phi = AR_coefs, sigma_sq = AR_resid_var, aic = AR_AIC)
}

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

#' Smooth AR coefficients and white noise variance
#'
#' @param spatial See \code{fit_bayesglm} internal code.
#' @param AR A Vxp matrix of estimated AR coefficients, where V is the number of vertices and p is the AR model order
#' @param var A vector length V containing the white noise variance estimates from the AR model
#' @param FWHM FWHM parameter for smoothing. Remember that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} or \code{NULL}
#'  to not do any smoothing. Default: \code{5}.
#'
#' @importFrom ciftiTools smooth_cifti make_surf as_xifti
#'
#' @keywords internal
#'
#' @return Smoothed AR coefficients and residual variance at every vertex
pw_smooth <- function(spatial, AR, var, FWHM=5){

  nV <- get_nV(spatial)
  if (nV$mdata != nrow(AR)) { stop('Number of rows in `AR` must match number of locations.') }
  if (nV$mdata != length(var)) { stop('Length of `var` must match number of locations.') }

  if (spatial$spatial_type == "vertex") {
    maskMdat <- if (is.null(spatial$maskMdat)) {
      rep(TRUE, nrow(spatial$surf$vertices))
    } else {
      spatial$maskMdat
    }

    # Just use left cortex, even though the data may be for the right cortex.

    AR_xif <- ciftiTools::as_xifti(
      cortexL = AR,
      surfL = spatial$surf,
      cortexL_mwall = maskMdat
    )
    #AR_xif$meta$cifti$brainstructures <- "left"
    AR_smoothed <- suppressWarnings(smooth_cifti(AR_xif, surf_FWHM = FWHM))
    AR_smoothed <- AR_smoothed$data$cortex_left

    var_xif <- ciftiTools::as_xifti(
      cortexL = var,
      surfL = spatial$surf,
      cortexL_mwall = maskMdat
    )
    #var_xif$meta$cifti$brainstructures <- "left"
    var_smoothed <- suppressWarnings(smooth_cifti(var_xif, surf_FWHM = FWHM))
    var_smoothed <- var_smoothed$data$cortex_left

  } else if (spatial$spatial_type == "voxel") {
    AR_xif <- ciftiTools::as_xifti(
      subcortVol = AR,
      subcortLabs = spatial$labels,
      subcortMask = spatial$maskMdat
    )
    AR_smoothed <- suppressWarnings(smooth_cifti(AR_xif, vol_FWHM = FWHM))
    AR_smoothed <- AR_smoothed$data$subcort

    var_xif <- ciftiTools::as_xifti(
      subcortVol = as.matrix(var),
      subcortLabs = spatial$labels,
      subcortMask = spatial$maskMdat
    )
    var_smoothed <- suppressWarnings(smooth_cifti(var_xif, vol_FWHM = FWHM))
    var_smoothed <- var_smoothed$data$subcort

  } else { stop() }

  list(AR = AR_smoothed, var = var_smoothed)
}
