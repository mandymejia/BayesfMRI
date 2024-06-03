#' Summarize a \code{"BGLM0"} object
#'
#' Summary method for class \code{"BGLM0"}
#'
#' @param object Object of class \code{"BGLM0"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BGLM0"} object, a list summarizing the properties
#'  of \code{object}.
#' @method summary BGLM0
summary.BGLM0 <- function(object, ...) {

  x <- list(
    fields = object$field_names,
    sessions = object$session_names,
    n_loc_total = length(object$spatial$mask),
    n_loc_modeled = sum(object$spatial$mask),
    GLM_type = attr(object$field_estimates, "GLM_type")
  )

  class(x) <- "summary.BGLM0"

  return(x)
}

#' @rdname summary.BGLM0
#' @export
#'
#' @param x Object of class \code{"summary.BGLM0"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BGLM0
print.summary.BGLM0 <- function(x, ...) {
  cat("====BGLM0 result====================\n")
  cat("Fields:   ", paste0("(", length(x$fields), ") ", paste(x$fields, collapse=", ")), "\n")
  if (length(x$sessions)==1 && x$sessions == "session_combined") {
    cat("Sessions: ", paste0("(", x$n_sess_orig, ", combined) \n"))
  } else {
    cat("Sessions: ", paste0("(", length(x$sessions), ") ", paste(x$sessions, collapse=", ")), "\n")
  }
  cat("Locations:", x$n_loc_modeled, "modeled,", x$n_loc_total, "total", "\n")
  cat("GLM type: ", x$GLM_type, "\n")
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.BGLM0
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BGLM0
print.BGLM0 <- function(x, ...) {
  print.summary.BGLM0(summary(x))
}
