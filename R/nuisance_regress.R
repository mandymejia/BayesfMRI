#' Regresses nuisance parameters from response or predictor variables
#'
#' @param Y a T by M matrix representing the response or predictor variables from which to regress nuisance.
#' @param Z a T by J matrix where J is the number of nuisance variables. It contains of nuisance regressors.
#'
#' @return A T by M matrix representing the residuals of Y after regressing out Z
#' @export
#'
#' @examples \dontrun{}
nuisance_regress <- function(Y, Z){

	if(nrow(Y) != nrow(Z)) stop('X and Z must have same number of rows')

 	invZtZ <- solve(t(Z) %*% Z)  #(Z'Z)^{-1}
	betahat <- invZtZ %*% t(Z) %*% Y #(Z'Z)^{-1} Z'Y
	Y2 <- Y - Z %*% betahat
	return(Y2)

}


