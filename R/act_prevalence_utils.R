#' Summarize a \code{"prev_BayesGLM"} object
#'
#' Summary method for class \code{"prev_BayesGLM"}
#'
#' @param object Object of class \code{"prev_BayesGLM"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.prev_BayesGLM"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary prev_BayesGLM
summary.prev_BayesGLM <- function(object, ...) {
  class(x) <- "summary.prev_BayesGLM"
  x
}

#' @rdname summary.prev_BayesGLM
#' @export
#'
#' @param x Object of class \code{"summary.prev_BayesGLM"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.prev_BayesGLM
print.summary.prev_BayesGLM <- function(x, ...) {
  cat("====BayesGLM Prevalences====================\n")
  cat("Summary for prevalences is not implemented yet.\n")
  invisible(NULL)
}

#' @rdname summary.prev_BayesGLM
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print prev_BayesGLM
print.prev_BayesGLM <- function(x, ...) {
  print.summary.prev_BayesGLM(summary(x))
}

#' Summarize a \code{"prev_BayesGLM_cifti"} object
#'
#' Summary method for class \code{"prev_BayesGLM_cifti"}
#'
#' @param object Object of class \code{"prev_BayesGLM_cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.prev_BayesGLM_cifti"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary prev_BayesGLM_cifti
summary.prev_BayesGLM_cifti <- function(object, ...) {
  class(x) <- "summary.prev_BayesGLM_cifti"
  x
}

#' @rdname summary.prev_BayesGLM_cifti
#' @export
#'
#' @param x Object of class \code{"summary.prev_BayesGLM_cifti"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.prev_BayesGLM_cifti
print.summary.prev_BayesGLM_cifti <- function(x, ...) {
  cat("====BayesGLM_cifti Prevalences==============\n")
  cat("Summary for prevalences is not implemented yet.\n")
  invisible(NULL)
}

#' @rdname summary.prev_BayesGLM_cifti
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print prev_BayesGLM_cifti
print.prev_BayesGLM_cifti <- function(x, ...) {
  print.summary.prev_BayesGLM_cifti(summary(x))
}
