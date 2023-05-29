#' Summarize a \code{"act_BayesGLM"} object
#'
#' Summary method for class \code{"act_BayesGLM"}
#'
#' @param object Object of class \code{"act_BayesGLM"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary act_BayesGLM
summary.act_BayesGLM <- function(object, ...) {
  act <- object$activations[!vapply(object$activations, is.null, 0)]
  if ("p_values" %in% names(act)) { act <- list(single_session=act) }
  x <- list(
    activations = lapply(act, function(x){apply(x$active, 2, table)}),
    method=object$method,
    alpha=object$alpha,
    threshold=object$threshold,
    correction=object$correction
    #excur_method
  )
  for (ii in seq(length(x$activations))) {
    colnames(x$activations[[ii]]) <- object$task_names
  }
  class(x) <- "summary.act_BayesGLM"
  x
}

#' @rdname summary.act_BayesGLM
#' @export
#'
#' @param x Object of class \code{"summary.act_BayesGLM"}.
#' @method print summary.act_BayesGLM
print.summary.act_BayesGLM <- function(x, ...) {
  x$activations <- x$activations$single_session
  cat("====BayesGLM Activations====================\n")
  cat(paste0("Activated locations (", sum(x$activations[,1]), " modeled locations):\n"))
  for (ii in seq(ncol(x$activations))) {
    cat(paste0(
      "    ", paste0(colnames(x$activations)[ii], ": ", x$activations["TRUE",ii]), "\n"
    ))
  }
  cat("GLM type:   ", x$method, "\n")
  cat("alpha:      ", x$alpha, "\n")
  cat("Threshold:  ", x$threshold, "\n")
  if (x$correction != "not applicable") {
    cat("Correction: ", x$correction, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.act_BayesGLM
#' @export
#'
#' @method print act_BayesGLM
print.act_BayesGLM <- function(x, ...) {
  print.summary.act_BayesGLM(summary(x))
}

#' Summarize a \code{"act_BayesGLM_cifti"} object
#'
#' Summary method for class \code{"act_BayesGLM_cifti"}
#'
#' @param object Object of class \code{"act_BayesGLM_cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary act_BayesGLM_cifti
summary.act_BayesGLM_cifti <- function(object, ...) {
  x <- summary.act_BayesGLM(object, ...)
  class(x) <- "summary.act_BayesGLM_cifti"
  x
}

#' @rdname summary.act_BayesGLM_cifti
#' @export
#'
#' @param x Object of class \code{"summary.act_BayesGLM_cifti"}.
#' @method print summary.act_BayesGLM_cifti
print.summary.act_BayesGLM_cifti <- function(x, ...) {
  cat("====BayesGLM_cifti Activations==============\n")
  cat(paste0("Activated locations:\n"))
  for (mm in seq(length(x$activations))) {
    cat(paste0(
      "  ", 
      names(x$activations)[mm], " (", 
      sum(x$activations[[mm]][,1]), 
      " modeled locations)\n"
    ))
    for (ii in seq(ncol(x$activations[[mm]]))) {
      cat(paste0(
        "    ", paste0(colnames(x$activations[[mm]])[ii], ": ", x$activations[[mm]]["TRUE",ii]), "\n"
      ))
    }
  }
  cat("GLM type:   ", x$method, "\n")
  cat("alpha:      ", x$alpha, "\n")
  cat("Threshold:  ", x$threshold, "\n")
  if (x$correction != "not applicable") {
    cat("Correction: ", x$correction, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.act_BayesGLM_cifti
#' @export
#'
#' @method print act_BayesGLM_cifti
print.act_BayesGLM_cifti <- function(x, ...) {
  print.summary.act_BayesGLM_cifti(summary(x))
}
