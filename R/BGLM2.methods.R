#' Summarize a \code{"fit_bglm2"} object
#'
#' Summary method for class \code{"fit_bglm2"}
#'
#' @param object Object of class \code{"fit_bglm2"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.fit_bglm2"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary fit_bglm2
summary.fit_bglm2 <- function(object, ...) {

  x <- list(
    n_contrasts = length(object$contrasts),
    fields = object$field_names,
    sessions = object$session_names,
    n_loc_total = vapply(lapply(object$model_results, '[[', "mask"), length, 0),
    n_loc_modeled = vapply(lapply(object$model_results, '[[', "mask"), sum, 0),
    excursion_type = object$excursion_type
  )
  class(x) <- "summary.fit_bglm2"
  return(x)
}

#' @rdname summary.fit_bglm2
#' @export
#'
#' @param x Object of class \code{"summary.fit_bglm2"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.fit_bglm2
print.summary.fit_bglm2 <- function(x, ...) {
  cat("====fit_bglm2 result===================\n")
  cat("Fields:   ", paste0("(", length(x$fields), ") ", paste(x$fields, collapse=", ")), "\n")
  if (length(x$sessions)==1 && x$sessions == "session_combined") {
    cat("Sessions: ", paste0("(", x$n_sess_orig, ", combined) \n"))
  } else {
    cat("Sessions: ", paste0("(", length(x$sessions), ") ", paste(x$sessions, collapse=", ")), "\n")
  }
  cat("Locations:\n")
  for (ii in seq(length(x$n_loc_total))) {
    cat(
      "          ", paste0(names(x$n_loc_total)[ii], ": ", x$n_loc_modeled[[ii]]),
      "modeled,", x$n_loc_total[[ii]], "total", "\n"
    )
  }
  cat("Excursion:", paste(x$excursion_type, collapse=", "), "\n")
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.fit_bglm2
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print fit_bglm2
print.fit_bglm2 <- function(x, ...) {
  print.summary.fit_bglm2(summary(x))
}

#' Summarize a \code{"BGLM2"} object
#'
#' Summary method for class \code{"BGLM2"}
#'
#' @param object Object of class \code{"BGLM2"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BGLM2"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary BGLM2
summary.BGLM2 <- function(object, ...) {
  x <- summary.fit_bglm2(object$fit_bglm2_results)
  class(x) <- "summary.BGLM2"
  x
}

#' @rdname summary.BGLM2
#' @export
#'
#' @param x Object of class \code{"summary.BGLM2"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BGLM2
print.summary.BGLM2 <- function(x, ...) {
  print.summary.fit_bglm2(x, ...)
}

#' @rdname summary.BGLM2
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BGLM2
print.BGLM2 <- function(x, ...) {
  print.summary.BGLM2(summary(x))
}
