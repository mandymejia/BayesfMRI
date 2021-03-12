#' Regresses nuisance parameters from response or predictor variables
#'
#' @param Y a T by M matrix representing the response or predictor variables from which to regress nuisance.
#' @param Z a T by J matrix where J is the number of nuisance variables. It contains of nuisance regressors.
#'
#' @return A T by M matrix representing the residuals of Y after regressing out Z
#' @export
#'
nuisance_regress <- function(Y, Z){
	# if(nrow(Y) != nrow(Z)) stop('X and Z must have same number of rows')
 	# invZtZ <- solve(t(Z) %*% Z)  #(Z'Z)^{-1}
	# betahat <- invZtZ %*% t(Z) %*% Y #(Z'Z)^{-1} Z'Y
	# Y2 <- Y - Z %*% betahat
	# return(Y2)

  # https://stackoverflow.com/questions/19100600/extract-maximal-set-of-independent-columns-from-a-matrix
  # https://stackoverflow.com/questions/39167204/in-r-how-does-one-extract-the-hat-projection-influence-matrix-or-values-from-an
  qrd <- qr(Z)
  Z <- Z[, qrd$pivot[seq_len(qrd$rank)], drop=FALSE]
  qrd <- qr(Z)
  Qd <- qr.Q(qrd)
  I_m_H <- diag(nrow(Z)) - (Qd %*% t(Qd))
  if (nrow(Y)==nrow(Z)) {
    return(I_m_H %*% Y)
  } else if (ncol(Y)==nrow(Z)) {
    return(Y %*% I_m_H)
  } else {
    stop("Y and Z are not of compatible dimensions.")
  }
}


