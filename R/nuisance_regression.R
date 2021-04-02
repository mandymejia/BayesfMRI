#' Hat matrix
#'
#' Get the hat matrix from a design matrix
#'
#' @param design The \eqn{T \times Q} design matrix
#'
#' @return The \eqn{T \times T} hat matrix
#' 
#' @keywords internal
hat_matrix <- function(design){
  design <- as.matrix(design)
  # https://stackoverflow.com/questions/19100600/extract-maximal-set-of-independent-columns-from-a-matrix
  # https://stackoverflow.com/questions/39167204/in-r-how-does-one-extract-the-hat-projection-influence-matrix-or-values-from-an
  qrd <- qr(design)
  design <- design[, qrd$pivot[seq_len(qrd$rank)], drop=FALSE]
  qrd <- qr(design)
  Qd <- qr.Q(qrd)
  Qd %*% t(Qd)
}

#' Nuisance regression
#'
#' Performs nuisance regression. The data and design matrix must both be
#'  centered, or an intercept must be included in the design matrix!
#'
#' @param Y The \eqn{T \times V} or \eqn{V \times T} data.
#' @param design The \eqn{T \times Q} matrix of nuisance regressors.
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
  I_m_H <- diag(nrow(design)) - hat_matrix(design)
  if (nrow(Y)==nrow(design)) {
    return(I_m_H %*% Y)
  } else if (ncol(Y)==nrow(design)) {
    return(Y %*% I_m_H)
  } else {
    stop("Y and design are not of compatible dimensions.")
  }
}