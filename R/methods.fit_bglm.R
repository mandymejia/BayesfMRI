#' Summarize a \code{"fit_bglm"} object
#'
#' Summary method for class \code{"fit_bglm"}
#'
#' @param object Object of class \code{"fit_bglm"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.fit_bglm"} object, a list summarizing the properties
#'  of \code{object}.
#' @method summary fit_bglm
summary.fit_bglm <- function(object, ...) {

  x <- list(
    fields = object$field_names,
    sessions = object$session_names,
    n_loc_In = sum(object$spatial$maskIn),
    n_loc_Mdat = sum(object$spatial$maskMdat),
    GLM_type = attr(object$field_estimates, "GLM_type")
  )

  class(x) <- "summary.fit_bglm"

  return(x)
}

#' @rdname summary.fit_bglm
#' @export
#'
#' @param x Object of class \code{"summary.fit_bglm"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.fit_bglm
print.summary.fit_bglm <- function(x, ...) {
  cat("====fit_bglm result====================\n")
  cat("Fields:   ", paste0("(", length(x$fields), ") ", paste(x$fields, collapse=", ")), "\n")
  if (length(x$sessions)==1 && x$sessions == "session_combined") {
    cat("Sessions: ", paste0("(", x$n_sess_orig, ", combined) \n"))
  } else {
    cat("Sessions: ", paste0("(", length(x$sessions), ") ", paste(x$sessions, collapse=", ")), "\n")
  }
  cat("Locations:", x$n_loc_Mdat, "modeled",
      paste0("(", x$n_loc_In, " input)\n"))
  cat("GLM type: ", x$GLM_type, "\n")
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.fit_bglm
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print fit_bglm
print.fit_bglm <- function(x, ...) {
  print.summary.fit_bglm(summary(x))
}
