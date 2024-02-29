#' Classical GLM
#'
#' Classical GLM for BayesGLM
#'
#' @param data,spatial,spatial_type,session_names,field_names,design_type See \code{BayesGLM}.
#' @param vcols_ss,nT_ss,nD,var_resid,sqrtInv_all See \code{BayesGLM}.
#' @param compute_SE compute SE?
#' @return A list of results
#' @keywords internal
GLM_classical <- function(
  BOLD, design, nV_D,
  session_names, field_names, design_type,
  vcols_ss, nT_ss, nD,
  var_resid, sqrtInv_all,
  do_pw, compute_SE=TRUE
  ){

  # OLD ~~~~~~~~~~
  # #setup for classical GLM
  # y_ss <- data_ss$BOLD #a vector (grouped by location)
  # X_ss <- do.call(cbind, data_ss$design) #cbind non-expanded design matrices for each field
  # valid_cols_bigX <- (colSums(abs(X_ss)) > 0) #because X_ss is a big sparse matrix, any missing fields will manifest as a block of zeros
  # valid_cols_bigX[is.na(valid_cols_bigX)] <- FALSE
  # X_ss <- X_ss[,valid_cols_bigX]

  # #perform classical GLM after any prewhitening
  # beta_hat_s <- SE_beta_hat_s <- matrix(NA, nV_all, nK)
  # XTX_inv <- try(Matrix::solve(Matrix::crossprod(X_ss)))
  # if (inherits(XTX_inv, "try-error")) {
  #   stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
  # }
  # coef_s <- as.matrix(XTX_inv %*% t(X_ss) %*% y_ss) #a vector of (estimates for location 1, estimates for location 2, ...)
  # coef_s_mat <- matrix(coef_s, nrow = nV, ncol = nK_ss) #re-form into a VxK matrix
  # beta_hat_s[mask==TRUE,cols_ss] <- coef_s_mat
  # resid_s <- t(matrix(y_ss - X_ss %*% coef_s, nrow = nT))
  # ~~~~~~~~~~~~~~

  y <- c(BOLD)
  X <- do.call(cbind, design[vcols_ss])
  nK_ss <- sum(vcols_ss) # without empty columns # [TO DO] integrate CompareGLM
  nK2 <- 0 # [TO DO] integrate CompareGLM

  XTX_inv <- try(Matrix::solve(Matrix::crossprod(X)))
  if (inherits(XTX_inv, "try-error")) {
    stop("There is some numerical instability in the design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
  }
  beta_hat <- XTX_inv %*% t(X) %*% y # coefs estimates for location 1, then 2, ...
  #beta_hat <- matrix(beta_hat, nrow = nV_D, ncol = nK_ss) #re-form into a VxK matrix
  resids <- t(matrix(y - c(X %*% beta_hat), nrow = nT_ss))

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
  } else {
    SE_beta_hat <- NULL
  }

  list(
    estimates = beta_hat,
    SE_estimates = SE_beta_hat,
    resids = resids,
    DOF = DOF_true
  )
}
