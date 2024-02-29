#' Summarize a \code{"BayesGLM_cifti"} object
#'
#' Summary method for class \code{"BayesGLM_cifti"}
#'
#' @param object Object of class \code{"BayesGLM_cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BayesGLM_cifti"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary BayesGLM_cifti
summary.BayesGLM_cifti <- function(object, ...) {

  x <- lapply(object$BayesGLM_results, summary)
  x <- x[!vapply(object$BayesGLM_results, is.null, FALSE)]
  x <- list(
    fields = x[[1]]$fields,
    sessions = x[[1]]$sessions,
    n_sess_orig = x[[1]]$n_sess_orig,
    n_loc_total = lapply(x, '[[', "n_loc_total"),
    n_loc_modeled = lapply(x, '[[', "n_loc_modeled"),
    #xii = summary(x$estimates_xii$classical[[1]]),
    GLM_type = x[[1]]$GLM_type,
    design_is_multiple = x[[1]]$design_is_multiple
  )
  class(x) <- "summary.BayesGLM_cifti"

  return(x)
}

#' @rdname summary.BayesGLM_cifti
#' @export
#'
#' @param x Object of class \code{"summary.BayesGLM_cifti"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BayesGLM_cifti
print.summary.BayesGLM_cifti <- function(x, ...) {
  cat("====BayesGLM_cifti result==============\n")
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
  if (x$design_is_multiple) {
    cat("GLM type: ", "classical (multiple designs)", "\n")
  } else {
    cat("GLM type: ", x$GLM_type, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.BayesGLM_cifti
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BayesGLM_cifti
print.BayesGLM_cifti <- function(x, ...) {
  print.summary.BayesGLM_cifti(summary(x))
}
