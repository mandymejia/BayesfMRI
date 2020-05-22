#' Dual Regression
#'
#' @param dat Subject-level fMRI data (TxV)
#' @param GICA Group-level independent components (QxV)
#' @param scale A logical value indicating whether the fMRI timeseries should be scaled by the image standard deviation (see Details).
#'
#' @return A list containing the subject-level independent components S (QxV), subject-level mixing matrix A (TxQ), and the row- and column- centered fMRI data (TxV)
#' @export
#' @importFrom matrixStats colVars
#'
dual_reg <- function(dat, GICA, scale=FALSE){

  ntime <- nrow(dat) #length of timeseries
  nvox <- ncol(dat) #number of data locations
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')
  if(nvox != ncol(GICA)) stop('The number of voxels in dat and GICA must match')

  Q <- nrow(GICA) #number of ICs
  if(Q > nvox) warning('More ICs than voxels. Are you sure?')
  if(Q > ntime) warning('More ICs than time points. Are you sure?')

  #center timeseries data across space and time and standardize scale
  dat_ctr <- scale_BOLD(dat, scale=scale)

  #center each group IC over voxels
  GICA_t <- scale(t(GICA), scale=FALSE)
  GICA <- t(GICA_t)

	#estimate A (IC timeseries)
	A <- dat_ctr %*% GICA_t %*% solve(GICA %*% GICA_t)
	#fix scale of timeseries (sd=1)
	sd_A <- sqrt(colVars(A))
	D <- diag(1/sd_A)
  A <- A %*% D

	#estimate S (IC maps)
	S <- solve(a=(t(A) %*% A), b=(t(A) %*% dat_ctr))

	#return result
	result <- list(S = S, A = A, dat_ctr = dat_ctr)
	return(result)

}
