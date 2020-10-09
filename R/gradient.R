#' Calculate gradient
#'
#' @param x  A vector or matrix. If a matrix, compute gradient of each column of x. 
#'
#' @return Gradient of x
#' 
#' @export
gradient <- function(x){

	isMatrix <- is.matrix(x)
	if(!isMatrix) x <- matrix(x, ncol=1)

	ncols <- ncol(x)
	dx <- x*0
	for(i in 1:ncols){

		xi <- x[,i]
		dx1 <- c(0, diff(xi))
		dx2 <- c(diff(xi),0)

		count <- c(1,rep(2, length(xi)-2),1) # Number to divide by when averaging dx1 and dx2
		dx[,i] <- (dx1 + dx2)/count
	}

	if(!isMatrix) dx <- as.vector(dx)
	return(dx)
}
