#' Summarize a \code{"BfMRI_design"} object
#'
#' Summary method for class \code{"BfMRI_design"}
#'
#' @param object Object of class \code{"BfMRI_design"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BfMRI_design"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary BfMRI_design
summary.BfMRI_design <- function(object, ...) {
  x <- object[c("field_names")]
  class(x) <- "summary.BfMRI_design"
  return(x)
}

#' @rdname summary.BfMRI_design
#' @export
#'
#' @param x Object of class \code{"summary.BfMRI_design"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BfMRI_design
print.summary.BfMRI_design <- function(x, ...) {
  cat("====BfMRI_design =======================\n")
  cat("Fields:", paste(x$field_names, collapse=", "), "\n")
  cat("----------------------------------------\n")
  invisible(NULL)
}

#' @rdname summary.BfMRI_design
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BfMRI_design
print.BfMRI_design <- function(x, ...) {
  print.summary.BfMRI_design(summary(x))
}