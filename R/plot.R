#' S3 method: use \code{\link[ciftiTools]{view_xifti_surface}} to plot a \code{"BayesGLM_cifti"} object
#'
#' @param x An object of class "BayesGLM_cifti"
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param method "Bayes" or "classical". \code{NULL} (default) will use
#'  the Bayesian results if available, and the classical results if not.
#' @param idx The data columns to plot. Overrides the \code{idx} argument to
#'  \code{\link[ciftiTools]{view_xifti_surface}}. \code{NULL} (default) will
#'  use all the data columns.
#' @param zlim Overrides the \code{zlim} argument for
#'  \code{\link[ciftiTools]{view_xifti_surface}}. Default: \code{c(-1, 1)}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti_surface}}
#'
#' @method plot BayesGLM_cifti
#'
# @importFrom ciftiTools view_xifti_surface
#' @export
#'
plot.BayesGLM_cifti <- function(x, session=NULL, method=NULL, idx=NULL, zlim=c(-1, 1), ...){

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("This function requires the `ciftiTools` package. Please install it.")
  }

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
  }
  if (is.null(x$task_estimates_xii[[method]][[session]])) {
    stop(paste("Session", session, "of method", method, "does not exist."))
  }

  # Column index
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, x$task_estimates_xii[[method]][[session]]$data)))
  }

  # Plot
  ciftiTools::view_xifti_surface(x$task_estimates_xii[[method]][[session]], idx=idx, zlim=zlim, ...)
}
