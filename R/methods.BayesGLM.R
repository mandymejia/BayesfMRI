#' Summarize a \code{"BGLM"} object
#'
#' Summary method for class \code{"BGLM"}
#'
#' @param object Object of class \code{"BGLM"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.BGLM"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary BGLM
summary.BGLM <- function(object, ...) {

  x <- lapply(object$BGLMs, summary)
  x <- x[!vapply(object$BGLMs, is.null, FALSE)]
  x <- list(
    fields = x[[1]]$fields,
    sessions = x[[1]]$sessions,
    n_sess_orig = x[[1]]$n_sess_orig,
    n_loc_Mdat = lapply(x, '[[', "n_loc_Mdat"),
    #xii = summary(x$estimate_xii$classical[[1]]),
    GLM_type = x[[1]]$GLM_type
  )
  class(x) <- "summary.BGLM"

  return(x)
}

#' @rdname summary.BGLM
#' @export
#'
#' @param x Object of class \code{"summary.BGLM"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.BGLM
print.summary.BGLM <- function(x, ...) {
  cat("====BGLM result==============\n")
  cat("Fields:   ", paste0("(", length(x$fields), ") ", paste(x$fields, collapse=", ")), "\n")
  if (length(x$sessions)==1 && x$sessions == "session_combined") {
    cat("Sessions: ", paste0("(", x$n_sess_orig, ", combined) \n"))
  } else {
    cat("Sessions: ", paste0("(", length(x$sessions), ") ", paste(x$sessions, collapse=", ")), "\n")
  }
  cat("Locations:\n")
  for (ii in seq(length(x$n_loc_Mdat))) {
    cat(
      "          ", paste0(names(x$n_loc_Mdat)[ii], ": ", x$n_loc_Mdat[[ii]]),
      "modeled data locations", "\n"
    )
  }
  cat("GLM type: ", x$GLM_type, "\n")
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.BGLM
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print BGLM
print.BGLM <- function(x, ...) {
  print.summary.BGLM(summary(x))
}
