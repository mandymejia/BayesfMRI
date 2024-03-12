#' Scale the BOLD timeseries
#'
#' @param BOLD fMRI data as a locations by time (\eqn{V \times T}) numeric 
#' 	matrix.
#' @param scale Option for scaling the BOLD response.
#' 
#' 	\code{"auto"} (default) will use \code{"mean"} scaling except if demeaned 
#'  data is detected (if any mean is less than one), in which case \code{"sd"}
#'  scaling will be used instead.
#' 
#' 	\code{"mean"} scaling will scale the data to percent local signal change.
#' 
#' 	\code{"sd"} scaling will scale the data by local standard deviation.
#' 
#' 	\code{"none"} will only center the data, not scale it. 
#'
#' @return Scale to units of percent local signal change and centers
#'
#' @importFrom stats var
#' @export
scale_BOLD <- function(BOLD, scale=c("auto", "mean", "sd", "none")){

	BOLD <- as.matrix(BOLD)
	scale <- match.arg(scale, c("auto", "mean", "sd", "none"))
	nV <- nrow(BOLD)
	nT <- ncol(BOLD)

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
			stop("Some local means are less than 0.1. Please set `scale_BOLD` to `'none'` or `'sd'`.")
		} else if (v_means_min < 1) {
			warning("Scaling to percent signal change when locations have means less than 1 may cause errors or produce aberrant results.")
		}
		BOLD <- 100*BOLD / v_means
	} else if (scale == "sd") {
		v_sd <- sqrt(apply(BOLD, 1, var, na.rm=TRUE))
		v_sd[is.na(v_sd)] <- 0
		if (min(v_sd) < 1e-6) {
			stop("Some local sds are less than 1e-6. Please set `scale_BOLD` to `'none'`.")
		}
		BOLD <- BOLD / v_sd
	}

	BOLD + v_means
}