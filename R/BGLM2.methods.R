#' Summarize a \code{"BGLM02"} object
#'
#' Summary method for class \code{"BGLM02"}
#'
#' @param object Object of class \code{"BGLM02"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BGLM02"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary BGLM02
summary.BGLM02 <- function(object, ...) {

  x <- list(
    n_contrasts = length(object$contrasts),
    fields = object$field_names,
    sessions = object$session_names,
    n_loc_total = vapply(lapply(object$model_results, '[[', "mask"), length, 0),
    n_loc_modeled = vapply(lapply(object$model_results, '[[', "mask"), sum, 0),
    excursion_type = object$excursion_type
  )
  class(x) <- "summary.BGLM02"
  return(x)
}

#' @rdname summary.BGLM02
#' @export
#'
#' @param x Object of class \code{"summary.BGLM02"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BGLM02
print.summary.BGLM02 <- function(x, ...) {
  cat("====BGLM02 result===================\n")
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

#' @rdname summary.BGLM02
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BGLM02
print.BGLM02 <- function(x, ...) {
  print.summary.BGLM02(summary(x))
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
  x <- summary.BGLM02(object$BGLM02_results)
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
  print.summary.BGLM02(x, ...)
}

#' @rdname summary.BGLM2
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BGLM2
print.BGLM2 <- function(x, ...) {
  print.summary.BGLM2(summary(x))
}
