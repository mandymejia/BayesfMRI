#' Estimate residual autocorrelation for prewhitening
#'
#' @param resids Estimated residuals
#' @param ar_order,aic Order of the AR model used to prewhiten the data at each location.
#'  If \code{!aic} (default), the order will be exactly \code{ar_order}. If \code{aic},
#'  the order will be between zero and \code{ar_order}, as determined by the AIC.
#' @importFrom stats ar.yw
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

#' Corrected AIC To-Do
#' 
#' Get corrected AIC
#' 
#' @keywords internal
cAIC <- function(...){invisible(NULL)}

#' Smooth AR coefficients and white noise variance
#'
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#' @param mask A logical vector indicating, for each vertex, whether to include
#'  it in smoothing. \code{NULL} (default) will use a vector of all \code{TRUE},
#'  meaning that no vertex is masked out; all are used for smoothing.
#' @param AR A Vxp matrix of estimated AR coefficients, where V is the number of vertices and p is the AR model order
#' @param var A vector length V containing the white noise variance estimates from the AR model
#' @param FWHM FWHM parameter for smoothing. Remember that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} or \code{NULL}
#'  to not do any smoothing. Default: \code{5}.#'
#'
#' @importFrom ciftiTools smooth_cifti make_surf
#'
#' @return Smoothed AR coefficients and residual variance at every vertex
pw_smooth <- function(vertices, faces, mask=NULL, AR, var, FWHM=5){

  if (is.null(mask)) { mask <- rep(TRUE, nrow(vertices)) }
  V <- sum(mask)
  V1 <- nrow(AR)
  V2 <- length(var)
  if(V != V1) stop('Number of rows in AR must match number of vertices')
  if(V != V2) stop('Length of var must match number of vertices')

  surf_smooth <- make_surf(
    list(
      pointset = vertices,
      triangle = faces
    )
  )
  AR_xif <- ciftiTools:::make_xifti(
    cortexL = AR,
    surfL = surf_smooth,
    cortexL_mwall = mask
  )
  #AR_xif$meta$cifti$brainstructures <- "left"
  AR_smoothed <- suppressWarnings(smooth_cifti(AR_xif, surf_FWHM = FWHM))
  AR_smoothed <- AR_smoothed$data$cortex_left

  var_xif <- ciftiTools:::make_xifti(
    cortexL = var,
    surfL = surf_smooth,
    cortexL_mwall = mask
  )
  #var_xif$meta$cifti$brainstructures <- "left"
  var_smoothed <- suppressWarnings(smooth_cifti(var_xif, surf_FWHM = FWHM))
  var_smoothed <- var_smoothed$data$cortex_left

  return(list(AR = AR_smoothed, var = var_smoothed))
}
