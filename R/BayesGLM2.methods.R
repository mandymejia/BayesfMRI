#' Summarize a \code{"BayesGLM2"} object
#'
#' Summary method for class \code{"BayesGLM2"}
#'
#' @param object Object of class \code{"BayesGLM2"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BayesGLM2"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary BayesGLM2
summary.BayesGLM2 <- function(object, ...) {

  x <- list(
    n_contrasts = length(object$contrasts),
    fields = object$field_names,
    sessions = object$session_names,
    n_loc_total = vapply(lapply(object$model_results, '[[', "mask"), length, 0),
    n_loc_modeled = vapply(lapply(object$model_results, '[[', "mask"), sum, 0),
    excursion_type = object$excursion_type
  )
  class(x) <- "summary.BayesGLM2"
  return(x)
}

#' @rdname summary.BayesGLM2
#' @export
#'
#' @param x Object of class \code{"summary.BayesGLM2"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BayesGLM2
print.summary.BayesGLM2 <- function(x, ...) {
  cat("====BayesGLM2 result===================\n")
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

#' @rdname summary.BayesGLM2
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BayesGLM2
print.BayesGLM2 <- function(x, ...) {
  print.summary.BayesGLM2(summary(x))
}

#' Summarize a \code{"BayesGLM2_cifti"} object
#'
#' Summary method for class \code{"BayesGLM2_cifti"}
#'
#' @param object Object of class \code{"BayesGLM2_cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BayesGLM2_cifti"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary BayesGLM2_cifti
summary.BayesGLM2_cifti <- function(object, ...) {
  x <- summary.BayesGLM2(object$BayesGLM2_results)
  class(x) <- "summary.BayesGLM2_cifti"
  x
}

#' @rdname summary.BayesGLM2_cifti
#' @export
#'
#' @param x Object of class \code{"summary.BayesGLM2_cifti"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BayesGLM2_cifti
print.summary.BayesGLM2_cifti <- function(x, ...) {
  print.summary.BayesGLM2(x, ...)
}

#' @rdname summary.BayesGLM2_cifti
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BayesGLM2_cifti
print.BayesGLM2_cifti <- function(x, ...) {
  print.summary.BayesGLM2_cifti(summary(x))
}
