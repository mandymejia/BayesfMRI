#' Scale the BOLD timeseries
#'
#' @param BOLD Input fMRI data (V x T)
#'
#' @return Scale to units of percent local signal change and centers
#' @export
#'
#' @examples \dontrun{}
scale_timeseries <- function(BOLD, scale=TRUE){

	BOLD <- as.matrix(BOLD)
	nvox <- nrow(BOLD)
	ntime <- ncol(BOLD)

	#check orientation, send warning message and transpose if necessary
	if(ntime > nvox){
		warning('More columns than rows. Transposing matrix so rows are data locations and columns are time points')
		BOLD <- t(BOLD)
		nvox <- nrow(BOLD)
		ntime <- ncol(BOLD)
	}

	local_means <- matrix(rowMeans(BOLD, na.rm=TRUE), nrow=nrow(BOLD), ncol=ncol(BOLD)) #the mean over time for each voxel (the mean image)
	if(scale) BOLD <- t(100*(BOLD - local_means)/local_means) #scale to units of pct local signal change AND CENTER
	if(!scale) BOLD <- t((BOLD - local_means)) #just center
	return(BOLD)

}
