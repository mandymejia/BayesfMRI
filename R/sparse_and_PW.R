#' Organize data for Bayesian GLM
#'
#' Transforms the usual TxV BOLD data matrix Y into vector form, and
#' 	the usual TxK design matrix X into big sparse matrix form for use in
#' 	Bayesian GLM.
#'
#' The Bayesian GLM requires \code{y} (a vector of length TV containing the BOLD data)
#' 	and \code{X_k} (a sparse TVxV matrix corresponding to the kth field regressor) for each field k.
#' 	The design matrices are combined as \code{A=cbind(X_1,...,X_K)}.
#'
#' @param BOLD,design,spatial,spde See \code{fit_bayesglm}.
#' @param field_names,design_type See \code{fit_bayesglm}.
#' @param valid_cols,nT,sqrtInv_all See \code{fit_bayesglm}.
#'
#' @return A list containing fields \code{y} and \code{A} (see Details)
#'
#' @details The Bayesian GLM requires \code{y} (a vector of length TV containing the BOLD data)
#' and \code{X_k} (a sparse TVxV matrix corresponding to the kth field regressor) for each field k.
#' The design matrices are combined as \code{A=cbind(X_1,...,X_K)}.
#'
#' @importFrom Matrix sparseMatrix
#'
#' @keywords internal
sparse_and_PW <- function(
  BOLD, design,
  spatial, spde,
  field_names, design_type,
  valid_cols, nT,
  sqrtInv_all
  ){
  nV <- get_nV(spatial)
  nK <- length(field_names)
	nIX <- seq(nT*nV$D)
	nIY <- rep(seq(nV$D), each = nT)

	y <- as.vector(BOLD) #makes a vector (y_1,...,y_V), where y_v is the timeseries for data location v
  X <- design

  A_sparse <- make_A_mat(spatial)

  ### Make `X_all` (design) and `bigX`. -----
	X_all <- vector('list', length=nK)
	for (kk in seq(nK)) {
	  if (design_type == "regular") {
	    # For missing fields.
	    if (!valid_cols[kk]) { # formerly if (is.nan(X[1,kk]))
	      X[,kk] <- rep(NA, length(X[,kk]))
	    }
	    # Expand the kth column of X into a VT x V.
	    # Then in `fit_bayesglm`: will post-multiply by A to get a VT x V2 matrix
	    #   (a V x V2 matrix for each time point).
	    X_all[[kk]] <- Matrix::sparseMatrix(nIX, nIY, x=rep(X[,kk], times=nV$D))
	    # #needs to be c-binded before model fitting.  For Bayesian GLM, post-multiply each by A before cbind().

	    # # previous approach
	    # X_k <- Matrix::sparseMatrix(nIX, nIY, x=rep(X[,kk], nV$D)) # %*% A #multiply by A to expand to the non-data locations
	    # bigX <- if (kk==1) { X_k } else { cbind(bigX, X_k) }
	  } else if (design_type == "per_location") {
	    if (!valid_cols[kk]) { # formerly if (is.nan(X[1,kk]))
	      X[,kk,] <- rep(NA, prod(dim(X)[c(1,3)]))
	    }
	    X_all[[kk]] <- Matrix::sparseMatrix(nIX, nIY, x=c(X[,kk,]))
	  } else { stop() }
	}

  # Prewhiten, if applicable. -----
  if (!is.null(sqrtInv_all)) {
    y <- as.vector(sqrtInv_all %*% y)
    X_all <- lapply(X_all, function(X_all_kk) { sqrtInv_all %*% X_all_kk } )
  }

  # Return results. -----
	list(BOLD=y, design=X_all, A_sparse=A_sparse)
}
