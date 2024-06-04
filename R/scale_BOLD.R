#' Scale the BOLD timeseries
#'
#' @param BOLD fMRI data as a locations by time (\eqn{V \times T}) numeric
#' 	matrix.
#' @param scale Option for scaling the BOLD response.
#' @param v_means Original means of the BOLD data. ONLY provide if data has already been centered.
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
scale_BOLD <- function(BOLD, scale=c("mean", "sd", "none"), v_means = NULL){

	BOLD <- as.matrix(BOLD)
	scale <- match.arg(scale, c("mean", "sd", "none"))
	nV <- nrow(BOLD)
	nT <- ncol(BOLD)

	# Get `v_means`, the mean over time for each location (the mean image)
	if(is.null(v_means)){
	  v_means <- rowMeans(BOLD, na.rm=TRUE)
	  BOLD <- BOLD - v_means
	} else {
	  if(length(v_means) != nV) stop('`v_means` must be of length nrow(BOLD)')
	  if(!is.numeric(v_means)) stop('`v_means` must be numeric')
	}
	v_means_min <- min(v_means, na.rm = TRUE)

	# # Determine `"auto"` scaling.
	# if (scale == "auto") {
	# 	scale <- if (v_means_min > 1) { "mean" } else { "sd" }
	# 	# cat("Using", scale, "scaling.\n")
	# }

	# Center and scale.
	if (scale == "mean") {
		if (v_means_min < .1) {
			stop("Some local means are less than 0.1. Please set `scale_BOLD` to `'none'` or `'sd'`,  or adjust `meanTol` to exclude low-mean locations.")
		} else if (v_means_min < 1) {
			warning("Scaling to percent signal change when locations have means less than 1 may cause errors or produce aberrant results. Consider setting `scale_BOLD` to `'none'` or `'sd'`,  or adjusting `meanTol` to exclude low-mean locations.")
		}
		BOLD <- 100* (BOLD / v_means) #since BOLD has already been centered, we end up with (BOLD - vmeans)/vmeans * 100 which represents % local signal change
	} else if (scale == "sd") {
		v_sd <- sqrt(apply(BOLD, 1, var, na.rm=TRUE))
		v_sd[is.na(v_sd)] <- 0
		if (min(v_sd) < 1e-6) {
			stop("Some local sds are less than 1e-6. Please set `scale_BOLD` to `'none'` or adjust `varTol` to exclude low-variance locations.")
		}
		BOLD <- BOLD / v_sd
	}

	return(BOLD)
}
