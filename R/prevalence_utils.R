#' Summarize a \code{"prev_fit_bglm"} object
#'
#' Summary method for class \code{"prev_fit_bglm"}
#'
#' @param object Object of class \code{"prev_fit_bglm"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.prev_fit_bglm"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary prev_fit_bglm
summary.prev_fit_bglm <- function(object, ...) {
  class(x) <- "summary.prev_fit_bglm"
  x
}

#' @rdname summary.prev_fit_bglm
#' @export
#'
#' @param x Object of class \code{"summary.prev_fit_bglm"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.prev_fit_bglm
print.summary.prev_fit_bglm <- function(x, ...) {
  cat("====BayesGLM Prevalences====================\n")
  cat("Summary for prevalences is not implemented yet.\n")
  invisible(NULL)
}

#' @rdname summary.prev_fit_bglm
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print prev_fit_bglm
print.prev_fit_bglm <- function(x, ...) {
  print.summary.prev_fit_bglm(summary(x))
}

#' Summarize a \code{"prev_BGLM"} object
#'
#' Summary method for class \code{"prev_BGLM"}
#'
#' @param object Object of class \code{"prev_BGLM"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.prev_BGLM"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary prev_BGLM
summary.prev_BGLM <- function(object, ...) {
  class(x) <- "summary.prev_BGLM"
  x
}

#' @rdname summary.prev_BGLM
#' @export
#'
#' @param x Object of class \code{"summary.prev_BGLM"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.prev_BGLM
print.summary.prev_BGLM <- function(x, ...) {
  cat("====BayesGLM Prevalences==============\n")
  cat("Summary for prevalences is not implemented yet.\n")
  invisible(NULL)
}

#' @rdname summary.prev_BGLM
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print prev_BGLM
print.prev_BGLM <- function(x, ...) {
  print.summary.prev_BGLM(summary(x))
}
