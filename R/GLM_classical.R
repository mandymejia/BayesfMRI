#' Classical GLM
#'
#' Classical GLM for BayesGLM
#'
#' @param BOLD,design,nK2,nV_D,field_names,design_type See \code{BayesGLM}.
#' @param vcols_ss,nT_ss,nD,var_resid,sqrtInv_all,do_pw See \code{BayesGLM}.
#' @param compute_SE compute SE?
#' @return A list of results
#' @keywords internal
GLM_classical <- function(
  BOLD, design, nK2, nV_D,
  field_names, design_type,
  vcols_ss, nT_ss, nD,
  var_resid, sqrtInv_all,
  do_pw, compute_SE=TRUE
  ){

  y <- c(BOLD)
  X <- do.call(cbind, design[vcols_ss])
  nK_ss <- sum(vcols_ss) # without empty columns # [TO DO] integrate CompareGLM

  XTX_inv <- try(Matrix::solve(Matrix::crossprod(X)))
  if (inherits(XTX_inv, "try-error")) {
    stop("There is some numerical instability in the design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
  }
  coefs <- XTX_inv %*% t(X) %*% y # coefs estimates for location 1, then 2, ...
  beta_hat <- matrix(coefs, nrow = nV_D, ncol = nK_ss) #re-form into a VxK matrix
  resids <- t(matrix(as.matrix(y) - (X %*% coefs), nrow = nT_ss))

  if (compute_SE) {
    # Residual SD.
    DOF_true <- (length(y)/nV_D) - nK_ss - nK2 - 1
    DOF_false <- (length(y)/nV_D - 1)
    var_error <- matrixStats::rowVars(resids) * DOF_false / DOF_true #correct DOF
    if(do_pw) var_error <- rep(mean(var_error), length(var_error)) #if prewhitening has been done, use same estimate of residual SD everywhere
    sd_error <- sqrt(var_error) #length = nV_D
    # SE of betas.
    # blocks of XTX_inv represent fields, so we should repeat each location-specific SD K times
    SE_beta_hat <- sqrt(Matrix::diag(XTX_inv)) * rep(sd_error, each = nK_ss)
    SE_beta_hat <- matrix(SE_beta_hat, ncol = nK_ss)
  } else {
    SE_beta_hat <- NULL
  }

  colnames(beta_hat) <- colnames(SE_beta_hat) <- field_names[vcols_ss]

  list(
    estimates = beta_hat,
    SE_estimates = SE_beta_hat,
    resids = resids,
    DOF = DOF_true,
    valid_fields = setNames(vcols_ss, field_names)
  )
}
