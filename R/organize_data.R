#' Organize data for Bayesian GLM
#'
#' Transforms the usual TxV BOLD data matrix Y into vector form, and
#' 	the usual TxK design matrix X into big sparse matrix form for use in
#' 	Bayesian GLM.
#'
#' The Bayesian GLM requires \code{y} (a vector of length TV containing the BOLD data)
#' 	and \code{X_k} (a sparse TVxV matrix corresponding to the kth task regressor) for each task k.
#' 	The design matrices are combined as \code{A=cbind(X_1,...,X_K)}.
#'
#' @param y the TxV data matrix containing the fMRI timeseries
#' @param X the TxK design matrix with K task-related columns
#' @param transpose Check orientation of data, which, if \code{TRUE}, will transpose
#' 	the data when the number of time points is greater than the number of voxels.
#' 	Note: this is not always true for subcortical regions.
#'
#' @return A list containing fields \code{y} and \code{A} (see Details)
#'
#' @details The Bayesian GLM requires \code{y} (a vector of length TV containing the BOLD data)
#' and \code{X_k} (a sparse TVxV matrix corresponding to the kth task regressor) for each task k.
#' The design matrices are combined as \code{A=cbind(X_1,...,X_K)}.
#'
#' @importFrom Matrix sparseMatrix
#'
#' @export
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

	result <- list(BOLD=y, design=bigX)
	return(result)
}

#' Organize prewhitened data for Bayesian GLM
#'
#' Transforms the usual TxV BOLD data matrix Y into vector form, and
#' 	the usual TxK design matrix X into big sparse matrix form for use in
#' 	Bayesian GLM.
#'
#' The Bayesian GLM requires \code{y} (a vector of length TV containing the BOLD data)
#' 	and \code{X_k} (a sparse TVxV matrix corresponding to the kth task regressor) for each task k.
#' 	The design matrices are combined as \code{A=cbind(X_1,...,X_K)}.
#'
#' @param y the TxV data matrix containing the fMRI timeseries
#' @param X the TxK design matrix with K task-related columns
#' @param transpose Check orientation of data, which, if \code{TRUE}, will transpose
#' 	the data when the number of time points is greater than the number of voxels.
#' 	Note: this is not always true for subcortical regions.
#'
#' @return A list containing fields \code{y} and \code{A} (see Details)
#'
#' @details The Bayesian GLM requires \code{y} (a vector of length TV containing the BOLD data)
#' and \code{X_k} (a sparse TVxV matrix corresponding to the kth task regressor) for each task k.
#' The design matrices are combined as \code{A=cbind(X_1,...,X_K)}.
#'
#' @importFrom Matrix bdiag
#'
#' @export
organize_data_pw <- function(y, X, transpose = TRUE){

  ntime <- nrow(y)
  nvox <- ncol(y)
  is_missing <- is.na(y[1,])
  K <- ncol(X) / sum(!is_missing)
  V <- sum(!is_missing)

  #check orientation, send warning message and transpose if necessary
  if(ntime > nvox & transpose == TRUE){
    warning('More columns than rows. Transposing matrix so rows are data locations and columns are time points')
    y <- t(y)
    ntime <- nrow(y)
    nvox <- ncol(y)
  }

  y <- as.vector(y) #makes a vector (y_1,...,y_V), where y_v is the timeseries for data location v

  # bigX <- Matrix(0,ntime*V,V*K)
  for(k in 1:K){
    k_inds <- seq(k,V*K,by = K)
    X_k <- X[,k_inds]
    if(k==1) bigX <- X_k else bigX <- cbind(bigX, X_k)
  }

  result <- list(y=y, X=bigX)
  return(result)
}
