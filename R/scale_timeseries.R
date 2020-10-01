#' Scale the BOLD timeseries
#'
#' @param BOLD Input fMRI data (V x T)
#' @param scale Option for scaling units to percent local signal change AND centering. If false, only centers.
#' @param transpose Check orientation of data, which, if TRUE, will transpose the data when the number of time points is greater than the number of voxels. Note: this is not always true for subcortical regions.
#'
#' @return Scale to units of percent local signal change and centers
#' @export
scale_timeseries <- function(BOLD, scale=TRUE, transpose = TRUE){

	BOLD <- as.matrix(BOLD)
	nvox <- nrow(BOLD)
	ntime <- ncol(BOLD)

	#check orientation, send warning message and transpose if necessary
	if(ntime > nvox & transpose == TRUE){
		warning('More columns than rows. Transposing matrix so rows are data locations and columns are time points')
		BOLD <- t(BOLD)
		nvox <- nrow(BOLD)
		ntime <- ncol(BOLD)
	}

	local_means <- matrix(rowMeans(BOLD, na.rm=TRUE), nrow=nrow(BOLD), ncol=ncol(BOLD)) #the mean over time for each voxel (the mean image)

	if(min(local_means, na.rm = T) < 1) warning("Scaling to percent signal change when locations have means less than 1 may cause errors or produce aberrant results.")
	if(min(local_means, na.rm = T) < 0.1) stop("Some local means are less than 0.1. Please set scale_BOLD = FALSE.")
	if(scale) BOLD <- t(100*(BOLD - local_means)/local_means) #scale to units of pct local signal change AND CENTER
	if(!scale) BOLD <- t((BOLD - local_means)) #just center
	return(BOLD)

}
