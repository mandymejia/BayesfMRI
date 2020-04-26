#' Transforms data and design matrices into the Bayesian GLM format
#'
#' @description Transforms the usual TxV BOLD data matrix Y into vector form, and
#' the usual TxK design matrix X into big sparse matrix form.
#'
#' @param y the TxV data matrix containing the fMRI timeseries
#' @param X the TxK design matrix with K task-related columns
#' @param transpose Check orientation of data, which, if TRUE, will transpose the data when the number of time points is greater than the number of voxels. Note: this is not always true for subcortical regions.
#'
#' @return A list containing fields `y` and `A` (see Details)
#'
#' @details The Bayesian GLM requires `y` (a vector of length TV containing the BOLD data)
#' and `X_k` (a sparse TVxV matrix corresponding to the kth task regressor) for each task k.
#' The design matrices are combined as `A=cbind(X_1,...,X_K)`.
#'
#' @export
#' @importFrom Matrix sparseMatrix
#'
#'
organize_data <- function(y, X, transpose = TRUE){

	ntime <- nrow(y)
	nvox <- ncol(y)

	#check orientation, send warning message and transpose if necessary
	if(ntime > nvox & transpose == TRUE){
		warning('More columns than rows. Transposing matrix so rows are data locations and columns are time points')
		y <- t(y)
		ntime <- nrow(y)
		nvox <- ncol(y)
	}

	y <- as.vector(y) #makes a vector (y_1,...,y_V), where y_v is the timeseries for data location v
	ix <- 1:(ntime*nvox)
	iy <- rep(1:nvox, each = ntime)

	K <- ncol(X)
	for(k in 1:K){
		X_k <- Matrix::sparseMatrix(ix, iy, x=rep(X[,k], nvox))
		if(k==1) bigX <- X_k else bigX <- cbind(bigX, X_k)
	}

	result <- list(y=y, X=bigX)
	return(result)
}
