#' Summarize a \code{"BayesGLM"} object
#'
#' Summary method for class \code{"BayesGLM"}
#'
#' @param object Object of class \code{"BayesGLM"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BayesGLM"} object, a list summarizing the properties
#'  of \code{object}.
#' @method summary BayesGLM
summary.BayesGLM <- function(object, ...) {

  x <- list(
    fields = object$field_names,
    sessions = object$session_names,
    n_sess_orig = object$n_sess_orig,
    n_loc_total = length(object$mask),
    n_loc_modeled = sum(object$mask),
    GLM_type = attr(object$field_estimates, "GLM_type"),
    design_is_multiple = !is.null(object$result_multiple)
  )

  class(x) <- "summary.BayesGLM"

  return(x)
}

#' @rdname summary.BayesGLM
#' @export
#'
#' @param x Object of class \code{"summary.BayesGLM"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BayesGLM
print.summary.BayesGLM <- function(x, ...) {
  cat("====BayesGLM result====================\n")
  cat("Fields:   ", paste0("(", length(x$fields), ") ", paste(x$fields, collapse=", ")), "\n")
  if (length(x$sessions)==1 && x$sessions == "session_combined") {
    cat("Sessions: ", paste0("(", x$n_sess_orig, ", combined) \n"))
  } else {
    cat("Sessions: ", paste0("(", length(x$sessions), ") ", paste(x$sessions, collapse=", ")), "\n")
  }
  cat("Locations:", x$n_loc_modeled, "modeled,", x$n_loc_total, "total", "\n")
  if (x$design_is_multiple) {
    cat("GLM type: ", "classical (multiple designs)", "\n")
  } else {
    cat("GLM type: ", x$GLM_type, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.BayesGLM
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BayesGLM
print.BayesGLM <- function(x, ...) {
  print.summary.BayesGLM(summary(x))
}
