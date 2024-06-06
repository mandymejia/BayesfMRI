#' Classical GLM
#'
#' Classical GLM for BayesGLM_fun (internal function)
#'
#' @param BOLD BOLD timeseries in vector form (TVx1), result of \code{sparse_and_PW}
#' @param design List of large sparse design matrices (each TVxV), one per regressor, result of \code{sparse_and_PW}
#' @param nK2,nV_D,field_names,design_type See \code{BayesGLM_fun}.
#' @param valid_cols,nT See \code{BayesGLM_fun}.
#' @param do_pw Has prewhitening been performed on the data and design?
#' @param compute_SE Compute SE of model coefficients?
#' @return A list of results
#' @keywords internal
GLM_classical <- function(
  BOLD, design, nK2, nV_D,
  field_names, design_type,
  valid_cols, nT,
  do_pw, compute_SE=TRUE
  ){

  y <- c(BOLD)
  X <- do.call(cbind, design[valid_cols]) #note that design is a list after running sparse_and_PW
  nK <- sum(valid_cols) # without empty columns # [TO DO] integrate CompareGLM

  XTX_inv <- try(Matrix::solve(Matrix::crossprod(X)))
  if (inherits(XTX_inv, "try-error")) {
    stop("There is some numerical instability in the design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
  }
  coefs <- XTX_inv %*% t(X) %*% y # coefs estimates for location 1, then 2, ...
  beta_hat <- matrix(coefs, nrow = nV_D, ncol = nK) #re-form into a VxK matrix
  resids <- t(matrix(as.matrix(y) - (X %*% coefs), nrow = nT))

  if (compute_SE) {
    # Residual SD.
    DOF_true <- (length(y)/nV_D) - nK - nK2 - 1
    DOF_false <- (length(y)/nV_D - 1)
    var_error <- matrixStats::rowVars(resids) * DOF_false / DOF_true #correct DOF
    if(do_pw) var_error <- rep(mean(var_error), length(var_error)) #if prewhitening has been done, use same estimate of residual SD everywhere
    sd_error <- sqrt(var_error) #length = nV_D
    # SE of betas.
    # blocks of XTX_inv represent fields, so we should repeat each location-specific SD K times
    SE_beta_hat <- sqrt(Matrix::diag(XTX_inv)) * rep(sd_error, each = nK)
    SE_beta_hat <- matrix(SE_beta_hat, ncol = nK)
    colnames(SE_beta_hat) <- field_names[valid_cols]
  } else {
    SE_beta_hat <- DOF_true <- NULL
  }

  colnames(beta_hat) <- field_names[valid_cols]

  list(
    estimates = beta_hat,
    SE_estimates = SE_beta_hat,
    resids = resids,
    RSS = rowSums(resids^2),
    DOF = DOF_true,
    valid_fields = setNames(valid_cols, field_names)
  )
}
