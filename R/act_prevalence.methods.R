#' Summarize a \code{"prev_BGLM0"} object
#'
#' Summary method for class \code{"prev_BGLM0"}
#'
#' @param object Object of class \code{"prev_BGLM0"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.prev_BGLM0"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary prev_BGLM0
summary.prev_BGLM0 <- function(object, ...) {
  class(x) <- "summary.prev_BGLM0"
  x
}

#' @rdname summary.prev_BGLM0
#' @export
#'
#' @param x Object of class \code{"summary.prev_BGLM0"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.prev_BGLM0
print.summary.prev_BGLM0 <- function(x, ...) {
  cat("====BayesGLM Prevalences====================\n")
  invisible(NULL)
}

#' @rdname summary.prev_BGLM0
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print prev_BGLM0
print.prev_BGLM0 <- function(x, ...) {
  print.summary.prev_BGLM0(summary(x))
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
  cat("====BayesGLM Prevalences====================\n")
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
