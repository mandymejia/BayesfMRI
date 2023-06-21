#' S3 method: use \code{\link[ciftiTools]{view_xifti_surface}} to plot a \code{"BayesGLM_cifti"} object
#'
#' @param x An object of class "BayesGLM_cifti"
#' @param idx Which task should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all tasks. This argument overrides 
#'  the \code{idx} argument to \code{\link[ciftiTools]{view_xifti_surface}}.
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param method "Bayes" or "classical". \code{NULL} (default) will use
#'  the Bayesian results if available, and the classical results if not.
#' @param zlim Overrides the \code{zlim} argument for
#'  \code{\link[ciftiTools]{view_xifti_surface}}. Default: \code{c(-1, 1)}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti_surface}}
#'
#' @method plot BayesGLM_cifti
#'
#' @importFrom ciftiTools view_xifti_surface
#' @export
#' 
#' @return Result of the call to \code{ciftiTools::view_cifti_surface}.
#'
plot.BayesGLM_cifti <- function(x, idx=NULL, session=NULL, method=NULL, zlim=c(-1, 1), ...){

  # Method
  if (is.null(method)) {
    method <- ifelse(
      is.null(x$task_estimates_xii$Bayes[[1]]),
      "classical", "Bayes"
    )
  }
  method <- match.arg(method, c("classical", "Bayes"))
  if (is.null(x$task_estimates_xii[[method]])) {
    stop(paste("Method", gsub("betas_", "", method, fixed=TRUE), "does not exist."))
  }

  # Session
  if (is.null(session)) {
    if (length(x$task_estimates_xii[[method]]) > 1) { message("Plotting the first session.") }
    session <- 1
  } else if (is.numeric(session)) {
    stopifnot(session %in% seq(length(x$session_names)))
  }
  the_xii <- x$task_estimates_xii[[method]][[session]]
  if (is.null(the_xii)) { stop(paste("Session", session, "does not exist.")) }

  # Column index
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, the_xii$data)))
  } else if (is.character(idx)) {
    idx <- match(idx, the_xii$meta$cifti$names)
  }

  # Plot
  ciftiTools::view_xifti_surface(the_xii, idx=idx, zlim=zlim, ...)
}

#' S3 method: use \code{\link[ciftiTools]{view_xifti_surface}} to plot a \code{"act_BayesGLM_cifti"} object
#'
#' @param x An object of class "act_BayesGLM_cifti"
#' @param idx Which task should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all tasks. This argument overrides 
#'  the \code{idx} argument to \code{\link[ciftiTools]{view_xifti_surface}}.
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti_surface}}
#'
#' @method plot act_BayesGLM_cifti
#'
#' @importFrom ciftiTools view_xifti_surface
#' @export
#' 
#' @return Result of the call to \code{ciftiTools::view_cifti_surface}.
#'
plot.act_BayesGLM_cifti <- function(x, idx=NULL, session=NULL, ...){

  # Session
  if (is.null(session)) {
    if (length(x$activations_xii) > 1) { message("Plotting the first session.") }
    session <- 1
  } else if (is.numeric(session)) {
    stopifnot(session %in% seq(length(x$activations_xii)))
  }
  the_xii <- x$activations_xii[[session]]
  if (is.null(the_xii)) { stop(paste("Session", session, "does not exist.")) }

  # Column index
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, the_xii$data)))
  } else if (is.character(idx)) {
    idx <- match(idx, the_xii$meta$cifti$names)
  }

  # Plot
  ciftiTools::view_xifti_surface(the_xii, idx=idx, ...)
}

#' S3 method: use \code{\link[ciftiTools]{view_xifti_surface}} to plot a \code{"BayesGLM2_cifti"} object
#'
#' @param x An object of class "BayesGLM2_cifti"
#' @param idx Which contrast should be plotted? Give the numeric index. 
#'  \code{NULL} (default) will show all contrasts. This argument overrides 
#'  the \code{idx} argument to \code{\link[ciftiTools]{view_xifti_surface}}.
#' @param what Estimates of the \code{"contrasts"} (default), or their 
#'  thresholded \code{"activations"}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti_surface}}
#'
#' @method plot BayesGLM2_cifti
#'
#' @importFrom ciftiTools view_xifti_surface
#' @export
#' 
#' @return Result of the call to \code{ciftiTools::view_cifti_surface}.
#'
plot.BayesGLM2_cifti <- function(x, idx=NULL, what=c("contrasts", "activations"), ...){
  what <- match.arg(what, c("contrasts", "activations"))
  what <- switch(what, contrasts="contrast_estimates_xii", activations="activations_xii")
  the_xii <- x[[what]]
  if (what=="activations_xii") {
    if (is.null(the_xii)) {
      stop("No activations in `'BayesGLM2'` object. Specify `excursion_type` in the `BayesGLM2` call and re-run.")
    }
    names(the_xii$meta$cifti$labels) <- paste0(
      names(the_xii$meta$cifti$labels), ", '",
      x$BayesGLM2_results$excursion_type, "'"
    )
  }
  
  # Column index
  if (is.null(idx)) { idx <- seq_len(ncol(do.call(rbind, the_xii$data))) }

  # Plot
  ciftiTools::view_xifti_surface(the_xii, idx=idx, ...)
}
