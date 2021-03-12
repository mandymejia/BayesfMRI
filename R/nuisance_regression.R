#' Nuisance regression
#'
#' Performs nuisance regression. The data and design matrix must both be
#'  centered, or an intercept must be included in the design matrix!
#'
#' @param Y The TxV or VxT data.
#' @param design The TxQ matrix of nuisance regressors.
#'
#' @return The data after nuisance regression.
#' 
#' @export
nuisance_regression <- function(Y, design){
  # Z <- design
	# if(nrow(Y) != nrow(Z)) stop('Y and Z must have same number of rows')
 	# invZtZ <- solve(t(Z) %*% Z)  #(Z'Z)^{-1}
	# betahat <- invZtZ %*% t(Z) %*% Y #(Z'Z)^{-1} Z'Y
	# return(Y - Z %*% betahat)

  Y <- as.matrix(Y); design <- as.matrix(design)
  # https://stackoverflow.com/questions/19100600/extract-maximal-set-of-independent-columns-from-a-matrix
  # https://stackoverflow.com/questions/39167204/in-r-how-does-one-extract-the-hat-projection-influence-matrix-or-values-from-an
  qrd <- qr(design)
  design <- design[, qrd$pivot[seq_len(qrd$rank)], drop=FALSE]
  qrd <- qr(design)
  Qd <- qr.Q(qrd)
  I_m_H <- diag(nrow(design)) - (Qd %*% t(Qd))
  if (nrow(Y)==nrow(design)) {
    return(I_m_H %*% Y)
  } else if (ncol(Y)==nrow(design)) {
    return(Y %*% I_m_H)
  } else {
    stop("Y and design are not of compatible dimensions.")
  }
}