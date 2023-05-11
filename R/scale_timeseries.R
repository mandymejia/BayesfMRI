#' Scale the BOLD timeseries
#'
#' @param BOLD Input fMRI data (V x T)
#' @param scale Option for scaling units.
#'
#' 	If \code{"auto"} (default), will use mean scaling except if demeaned data
#' 	is detected, in which case sd scaling will be used instead.
#'
#' 	\code{"mean"} scaling will scale the data to percent local signal change.
#'
#' 	\code{"sd"} scaling will scale the data by local standard deviation.
#'
#' 	\code{"none"} will only center the data, not scale it.
#' @param transpose Check orientation of data, which, if TRUE, will transpose
#' 	the data when the number of time points is greater than the number of voxels.
#' 	Note: this is not always true for subcortical regions.
#'
#' @return Scale to units of percent local signal change and centers
#'
#' @examples
#' nT <- 30
#' nV <- 400
#' X <- matrix(rnorm(nV*nT), nrow=nV)
#' scale_timeseries(X)
#' 
#' @export
scale_timeseries <- function(BOLD, scale=c("auto", "mean", "sd", "none"), transpose = TRUE){

	BOLD <- as.matrix(BOLD)
	nvox <- nrow(BOLD)
	ntime <- ncol(BOLD)

	scale <- match.arg(scale, c("auto", "mean", "sd", "none"))

	# Check orientation, send warning message and transpose if necessary.
	if((ntime > nvox) & transpose){
		warning('More columns than rows. Transposing matrix so rows are data locations and columns are time points')
		BOLD <- t(BOLD)
		nvox <- nrow(BOLD)
		ntime <- ncol(BOLD)
	}

	# Get `v_means`, the mean over time for each location (the mean image)
	v_means <- rowMeans(BOLD, na.rm=TRUE)
	v_means_min <- min(v_means, na.rm = TRUE)

	# Determine `"auto"` scaling.
	if (scale == "auto") {
		scale <- if (v_means_min > 1) { "mean" } else { "sd" }
		# cat("Using", scale, "scaling.\n")
	}

	# Center and scale.
	BOLD <- BOLD - v_means
	if (scale == "mean") {
		if (v_means_min < .1) {
			stop("Some local means are less than 0.1. Please set `scale_BOLD = 'none'` or set `meanTol` > .1 in the call to BayesGLM.")
		} else if (v_means_min < 1) {
			warning("Scaling to percent signal change when locations have means less than 1 may cause errors or produce aberrant results.")
		}
		BOLD <- t(100*BOLD / v_means)
	} else if (scale == "sd") {
		v_sd <- sqrt(rowVars(BOLD, na.rm=TRUE))
		v_sd[is.na(v_sd)] <- 0
		if (min(v_sd) < 1e-6) {
			stop("Some local sds are less than 1e-6. Please set `scale_BOLD = 'none' or set `varTol` > 1e-3 in the call to BayesGLM.")
		}
		BOLD <- t(BOLD / v_sd)
	} else {
		BOLD <- t(BOLD)
	}

	BOLD
}
