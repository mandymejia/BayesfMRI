#' PCA-based Dimension Reduction and Prewhitening for ICA
#'
#' @description Performs dimension reduction and prewhitening based on probablistic PCA using SVD. If dimensionality is not specified, it is estimated using the method described in Minka (2008).
#'
#' @param X VxT fMRI timeseries data matrix, centered by columns and rows (columns are actually all that matter, but MATLAB implementation of Minka method also assumes rows have been centered (implicit in use of cov function))
#' @param Q Number of latent dimensions to estimate. If not specified, estimated using PESEL (Sobczyka et al. 2020).
#' @param Q_max Maximal number of principal components for automatic dimensionality selection with PESEL
#'
#' @return A list containing the dimension-reduced data (data_reduced, a VxQ matrix), prewhitening/dimension reduction matrix (H, a QxT matrix) and its (pseudo-)inverse (Hinv, a TxQ matrix), the noise variance (sigma_sq), the correlation matrix of the dimension-reduced data (C_diag, a QxQ matrix), and the dimensionality (Q)
#' @export
#' @import pesel
#'
dim_reduce = function(X, Q=NULL, Q_max=100){

  nvox = nrow(X) #number of brain locations
  ntime = ncol(X) #number of fMRI volumes (reduce this)
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')

  # check that X has been centered both ways
  tol = 1e-12
  if(max(colMeans(X)^2) > tol) stop('Columns of X must be centered')
  if(max(rowMeans(X)^2) > tol) warning('Rows of X should be centered')

  #determine PCA dimensionality
  if(is.null(Q)){
    pesel_X <- pesel(X, npc.max=Q_max, method='homogenous')
    Q <- pesel_X$nPCs
  }

  #perform dimension reduction
  XXt = t(X) %*% X / nvox
  svd_XXt = svd(XXt, nu=Q, nv=0)
  U = svd_XXt$u
  D1 = svd_XXt$d[1:Q]
  D2 = svd_XXt$d[(Q+1):length(svd_XXt$d)]

  #residual variance
  sigma_sq = mean(D2)

  #prewhitening matrix
  H = diag(1/sqrt(D1 - sigma_sq)) %*% t(U)
  H_inv = U %*% diag(sqrt(D1 - sigma_sq))

  #for residual variance after prewhitening
  C_diag = diag(H %*% t(H))

  #prewhitened data (transposed)
  X_new <- X %*% t(H)

  result = list(data_reduced=X_new, H=H, H_inv=H_inv, sigma_sq=sigma_sq, C_diag=C_diag, Q=Q)
  return(result)

}
