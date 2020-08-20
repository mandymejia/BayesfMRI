#' Dual Regression
#'
#' @param dat Subject-level fMRI data (VxT)
#' @param GICA Group-level independent components (VxQ)
#' @param scale A logical value indicating whether the fMRI timeseries should be scaled by the image standard deviation.
#'
#' @return A list containing the subject-level independent components S (VxQ), subject-level mixing matrix A (TxQ), and the row- and column- centered fMRI data (VxT)
#' @export
#' @importFrom matrixStats colVars
#'
dual_reg <- function(dat, GICA, scale=FALSE){

  ntime <- ncol(dat) #length of timeseries
  nvox <- nrow(dat) #number of data locations
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')
  if(nvox != nrow(GICA)) stop('The number of voxels in dat and GICA must match')

  Q <- ncol(GICA) #number of ICs
  if(Q > nvox) warning('More ICs than voxels. Are you sure?')
  if(Q > ntime) warning('More ICs than time points. Are you sure?')

  #center timeseries data across space and time (and standardize scale if scale=TRUE)
  dat_ctr <- scale_BOLD(dat, scale=scale)
  dat_ctr_t <- t(dat_ctr)

  #center each group IC over voxels
  GICA <- scale(GICA, scale=FALSE)

	#estimate A (IC timeseries)
	A <- dat_ctr_t %*% GICA %*% solve(t(GICA) %*% GICA)
	#fix scale of timeseries (sd=1)
	sd_A <- sqrt(colVars(A))
	D <- diag(1/sd_A)
  A <- A %*% D

	#estimate S (IC maps)
	S <- solve(a=(t(A) %*% A), b=(t(A) %*% dat_ctr_t))

	#return result
	result <- list(S = S, A = A, dat_ctr = dat_ctr)
	return(result)

}
