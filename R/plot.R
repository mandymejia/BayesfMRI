#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"BGLM"} object
#'
#' @param x An object of class "BGLM"
#' @param Bayes \code{TRUE} for plotting Bayesian results, \code{FALSE} for plotting
#' classical GLM results. Default: \code{NULL}, which will use the Bayesian results
#' if available and the classical results if not.
#' @param idx Which field should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all fields. This argument overrides
#'  the \code{idx} argument to \code{\link[ciftiTools]{view_xifti}}.
#' @param title If NULL, the field names associated with idx will be used.
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param zlim Overrides the \code{zlim} argument for
#'  \code{\link[ciftiTools]{view_xifti}}. Default: \code{c(-1, 1)}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti}}
#'
#' @method plot BGLM
#'
#' @importFrom ciftiTools view_xifti
#' @export
#'
#' @return Result of the call to \code{ciftiTools::view_cifti}.
#'
plot.BGLM <- function(x, Bayes=NULL, idx=NULL, title=NULL, session=NULL, zlim=c(-1, 1), ...){

  # Method
  if (is.null(Bayes)) {
    method <- if (!is.null(x$estimate_xii$Bayes[[1]])) { "Bayes" } else { "classical" }
  } else if (isTRUE(Bayes) || isFALSE(Bayes)) {
    method <- if (isTRUE(Bayes)) { "Bayes" } else { "classical" }
  } else {
    stop('`Bayes` must be `TRUE`, `FALSE` or `NULL`.')
  }
  # if (is.null(x$estimate_xii[[method]])) {
  #   stop(paste("Method", gsub("betas_", "", method, fixed=TRUE), "does not exist."))
  # }

  # Session
  if (is.null(session)) {
    if (length(x$estimate_xii[[method]]) > 1) { message("Plotting the first session.") }
    session <- 1
  } else if (is.numeric(session)) {
    stopifnot(length(session)==1)
    stopifnot(session %in% seq(length(x$session_names)))
  }
  the_xii <- x$estimate_xii[[method]][[session]]
  if (is.null(the_xii)) { stop(paste("Session", session, "does not exist.")) }

  # Column index
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, the_xii$data)))
  } else if (is.character(idx)) {
    idx <- match(idx, the_xii$meta$cifti$names)
  }

  # Names
  idx_names <- the_xii$meta$cifti$names[idx]

  # Title(s)
  if(is.null(title)){
    title <- idx_names
  }

  # Plot
  ciftiTools::view_xifti(the_xii, idx=idx, title=title, zlim=zlim, fname_suffix = idx_names, ...)
}

#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"act_BGLM"} object
#'
#' @param x An object of class "act_BGLM"
#' @param idx Which field should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all fields. This argument overrides
#'  the \code{idx} argument to \code{\link[ciftiTools]{view_xifti}}.
#' @param title If NULL, the field names associated with idx will be used.
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti}}
#'
#' @method plot act_BGLM
#'
#' @importFrom ciftiTools view_xifti
#' @export
#'
#' @return Result of the call to \code{ciftiTools::view_cifti_surface}.
#'
plot.act_BGLM <- function(x, idx=NULL, title=NULL, session=NULL, ...){

  # Session
  if (is.null(session)) {
    if (length(x$activations_xii) > 1) { message("Plotting the first session.") }
    session <- 1
  } else if (is.numeric(session)) {
    stopifnot(length(session)==1)
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

  # Names
  idx_names <- the_xii$meta$cifti$names[idx]

  # Title(s)
  if(is.null(title)){
    title <- idx_names
  }

  # Values and colors
  vals <- unique(c(as.matrix(the_xii)))
  vals <- setdiff(vals, 0)
  #for a single level (e.g., only a single gamma value was given to activations()), use red to color activations
  if(length(vals) == 1){
    the_xii <- convert_to_dlabel(the_xii, colors = "red",
                                 labels = rownames(the_xii$meta$cifti$labels[[1]]))
  }

  # Plot
  ciftiTools::view_xifti(the_xii, idx=idx, title=title, fname_suffix = idx_names, ...)
}

#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"BGLM2"} object
#'
#' @param x An object of class "BGLM2"
#' @param idx Which contrast should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all contrasts. This argument
#'  overrides the \code{idx} argument to \code{\link[ciftiTools]{view_xifti}}.
#' @param what Estimates of the \code{"contrasts"} (default), or their
#'  thresholded \code{"activations"}.
#' @param zlim Overrides the \code{zlim} argument for
#'  \code{\link[ciftiTools]{view_xifti}}. Default: \code{c(-1, 1)}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti}}
#'
#' @method plot BGLM2
#'
#' @importFrom ciftiTools view_xifti
#' @export
#'
#' @return Result of the call to \code{ciftiTools::view_cifti}.
#'
plot.BGLM2 <- function(x, idx=NULL, what=c("contrasts", "activations"), zlim=c(-1, 1), ...){
  what <- match.arg(what, c("contrasts", "activations"))
  what <- switch(what, contrasts="contrast_estimate_xii", activations="activations_xii")
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
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, the_xii$data)))
  } else if (is.character(idx)) {
    idx <- match(idx, the_xii$meta$cifti$names)
  }

  # Plot
  ciftiTools::view_xifti(the_xii, idx=idx, zlim=zlim, ...)
}

#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"prev_BGLM"} object
#'
#' @param x An object of class "prev_BGLM"
#' @param idx Which task should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all tasks. This argument overrides
#'  the \code{idx} argument to \code{\link[ciftiTools]{view_xifti}}.
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param drop_zeros Color locations without any activation across all results
#'  (zero prevalence) the same color as the medial wall? Default: \code{NULL} to
#'  drop the zeros if only one \code{idx} is being plotted.
#' @param colors,zlim See \code{\link[ciftiTools]{view_xifti}}. Here, the defaults
#'  are overrided to use the Viridis \code{"plasma"} color scale between
#'  \code{1/nA} and 1, where \code{nA} is the number of results in \code{x}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti}}
#'
#' @method plot prev_BGLM
#'
#' @importFrom ciftiTools view_xifti
#' @importFrom fMRItools is_1
#' @export
#'
#' @return Result of the call to \code{ciftiTools::view_cifti_surface}.
#'
plot.prev_BGLM <- function(
  x, idx=NULL, session=NULL,
  drop_zeros=NULL, colors="plasma",
  zlim=c(round(1/x$n_results-.005, 2), 1), ...){

  # Session
  if (is.null(session)) {
    if (length(x$prev_xii) > 1) { message("Plotting the first session.") }
    session <- 1
  } else if (is.numeric(session)) {
    stopifnot(length(session)==1)
    stopifnot(session %in% seq(length(x$prev_xii)))
  }
  the_xii <- x$prev_xii[[session]]
  if (is.null(the_xii)) { stop(paste("Session", session, "does not exist.")) }

  # Column index
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, the_xii$data)))
  } else if (is.character(idx)) {
    idx <- match(idx, the_xii$meta$cifti$names)
  }

  if (is.null(drop_zeros)) { drop_zeros <- length(idx) == 1 }
  stopifnot(is_1(drop_zeros, "logical"))
  if (drop_zeros) {
    the_xii <- move_to_mwall(the_xii, 0)
  }

  # Plot
  ciftiTools::view_xifti(the_xii, idx=idx, colors=colors, zlim=zlim, ...)
}
