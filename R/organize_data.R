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
organize_data <- function(y, X, n_mesh, inds, sqrtInv_all=NULL){

  # Check inputs. -----
  stopifnot((is.matrix(y) | is.data.frame(y)) && is.numeric(y))
  stopifnot((is.matrix(X) | is.data.frame(X)) && is.numeric(X))

  ### Transpose `y` if necessary. -----
  if (ncol(y) == nrow(X)) {
    warning('Transposing fMRI data (`y`) so rows are time points and columns are locations.')
    y <- t(y)
  }

  # Define constants. -----
	nT <- nrow(y)
	nV <- ncol(y) #the number of data locations
	nV2 <- n_mesh # the number of mesh locations, which may be a superset
	nIX <- seq(nT*nV)
	nIY <- rep(seq(nV), each = nT)
	nK <- ncol(X)

  # Make outputs. -----
  ### `y.` -----
	y <- as.vector(y) #makes a vector (y_1,...,y_V), where y_v is the timeseries for data location v

  ### `A`. -----
	A <- Matrix::Diagonal(nV2)[inds,]

  ### Make `X_all` (design) and `bigX`. -----
	X_all <- vector('list', length=nK)
	for (kk in seq(nK)) {
	  # expand the kth column of X into a VT x V, then post-multiply by A to get a VT x V2 matrix (a V x V2 matrix for each time point)
    if(is.nan(X[1,kk])) X[,kk] <- rep(NA, length(X[,kk])) #for tasks that are missing for a particular session
    #[TO DO] If X varies spatially, replace rep() below with c() to vectorize the spatially-varying design.  new line: X_all_k <- ifelse(..., rep, c)
	  X_all[[kk]] <- Matrix::sparseMatrix(nIX, nIY, x=rep(X[,kk], times=nV)) #needs to be c-binded before model fitting.  For Bayesian GLM, post-multiply each by A before cbind().

	  # # previous approach
	  # X_k <- Matrix::sparseMatrix(nIX, nIY, x=rep(X[,kk], nV)) # %*% A #multiply by A to expand to the non-data locations
    # bigX <- if (kk==1) { X_k } else { cbind(bigX, X_k) }
	}


  # Prewhiten, if applicable. -----
  if (!is.null(sqrtInv_all)) {
     y <- as.vector(sqrtInv_all %*% y)
     X_all <- lapply(X_all, function(X_all_kk) { sqrtInv_all %*% X_all_kk } )
  #   #bigX <- sqrtInv_all %*% bigX
  }

  # Return results. -----
	list(
	  BOLD=y, design=X_all,
	  #bigX = bigX,
	  A=A
	 )
}
