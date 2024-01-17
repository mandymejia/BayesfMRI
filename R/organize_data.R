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
#' @param n_mesh the number of mesh locations V2 >= V
#' @param inds the indices in the mesh that correspond to the V data locations
#'
#' @return A list containing fields \code{y} and \code{A} (see Details)
#'
#' @details The Bayesian GLM requires \code{y} (a vector of length TV containing the BOLD data)
#' and \code{X_k} (a sparse TVxV matrix corresponding to the kth task regressor) for each task k.
#' The design matrices are combined as \code{A=cbind(X_1,...,X_K)}.
#'
#' @importFrom Matrix sparseMatrix
#'
#' @keywords internal
organize_data <- function(y, X, n_mesh, inds){

  stopifnot((is.matrix(y) | is.data.frame(y)) && is.numeric(y))
  stopifnot((is.matrix(X) | is.data.frame(X)) && is.numeric(X))

  if (ncol(y) == nrow(X)) {
    warning('Transposing fMRI data (`y`) so rows are time points and columns are locations.')
    y <- t(y)
  }
	nT <- nrow(y)
	nV <- ncol(y) #the number of data locations
	nV2 <- n_mesh # the number of mesh locations, which may be a superset

	y <- as.vector(y) #makes a vector (y_1,...,y_V), where y_v is the timeseries for data location v
	ix <- seq(nT*nV)
	iy <- rep(seq(nV), each = nT)

	A <- Matrix::Diagonal(nV2)[inds,]

	nK <- ncol(X)
	X_all <- vector('list', length=nK)
	for (kk in seq(nK)) {
	  # expand the kth column of X into a VT x V, then post-multiply by A to get a VT x V2 matrix (a V x V2 matrix for each time point)
    if(is.nan(X[1,kk])) X[,kk] <- rep(NA, length(X[,kk]))
	  X_all[[kk]] <- Matrix::sparseMatrix(ix, iy, x=rep(X[,kk], nV)) #needs to be c-binded before model fitting.  For Bayesian GLM, post-multiply each by A before cbind().

	  #previous approach
	  X_k <- Matrix::sparseMatrix(ix, iy, x=rep(X[,kk], nV)) # %*% A #multiply by A to expand to the non-data locations
    bigX <- if (kk==1) { X_k } else { cbind(bigX, X_k) }
	}

	list(BOLD=y, design=X_all, bigX = bigX, A=A)
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
#' @keywords internal
organize_data_pw <- function(y, X, transpose = TRUE){

  if (nrow(y) > ncol(y) && transpose) {
    warning('More columns than rows. Transposing matrix so rows are data locations and columns are time points')
    y <- t(y)
  }
	nT <- nrow(y)

  is_missing <- is.na(y[1,])
  nK <- ncol(X) / sum(!is_missing)
  nV <- sum(!is_missing)

  y <- as.vector(y) #makes a vector (y_1,...,y_V), where y_v is the timeseries for data location v

  # bigX <- Matrix(0,nT*V,V*K)
  for (kk in seq(nK)) {
    k_inds <- seq(kk, nV*nK,by = nK)
    X_k <- X[,k_inds]
    bigX <- if (kk==1) { X_k } else { cbind(bigX, X_k) }
  }

  list(y=y, X=bigX)
}
