#Z is the matrix of nuisance regressors
#Y is the response or predictor variables from which to regress nuisance

nuisance_regress <- function(Y, Z){

	if(nrow(Y) != nrow(Z)) stop('X and Z must have same number of rows')

 	invZtZ <- solve(t(Z) %*% Z)  #(Z'Z)^{-1}
	betahat <- invZtZ %*% t(Z) %*% Y #(Z'Z)^{-1} Z'Y
	Y2 <- Y - Z %*% betahat
	return(Y2)

}


