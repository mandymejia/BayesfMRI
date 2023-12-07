#' Summarize a \code{"act_BayesGLM"} object
#'
#' Summary method for class \code{"act_BayesGLM"}
#'
#' @param object Object of class \code{"act_BayesGLM"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.act_BayesGLM"} object, a list summarizing the 
#'  properties of \code{object}.
#' @method summary act_BayesGLM
summary.act_BayesGLM <- function(object, ...) {
  act <- object$activations[!vapply(object$activations, is.null, 0)]
  if ("p_values" %in% names(act)) { act <- list(single_session=act) }
  x <- list(
    activations = lapply(act, function(x_a){
      if ("active" %in% names(x_a)) {
        q <- apply(x_a$active, 2, function(avec){c(`TRUE`=sum(avec), `FALSE`=sum(!avec))})
        colnames(q) <- object$task_names
        q
      } else {
        lapply(x_a, function(x_b){
          q <- apply(x_b$active, 2, function(avec){c(`TRUE`=sum(avec), `FALSE`=sum(!avec))})
          colnames(q) <- object$task_names
          q
        })
      }
    }),
    method=object$method,
    alpha=object$alpha,
    gamma=object$gamma,
    correction=object$correction
    #excur_method
  )
  class(x) <- "summary.act_BayesGLM"
  x
}

#' @rdname summary.act_BayesGLM
#' @export
#'
#' @param x Object of class \code{"summary.act_BayesGLM"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.act_BayesGLM
print.summary.act_BayesGLM <- function(x, ...) {
  cat("====BayesGLM Activations====================\n")
  cat(paste0("Activated locations (", sum(x$activations[[1]][,1]), " modeled locations):\n"))
  for (ii in seq(length(x$activations))) {
    act_ii <- x$activations[[ii]]
    cat(paste0(
      "  ", names(x$activations)[ii],
      "\n"
    ))
    for (kk in seq(ncol(act_ii))) {
      cat(paste0(
        "    ", colnames(act_ii)[kk],
        ": ", act_ii["TRUE",kk], "\n"
      ))
    }
  }
  cat("GLM type:   ", x$method, "\n")
  cat("alpha:      ", x$alpha, "\n")
  cat("gamma:      ", x$gamma, "\n")
  if (x$correction != "not applicable") {
    cat("Correction: ", x$correction, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.act_BayesGLM
#' @export
#'
#' @return \code{NULL}, invisibly.
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
#' @return A \code{"summary.act_BayesGLM_cifti"} object, a list summarizing the 
#'  properties of \code{object}.
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
#' @return \code{NULL}, invisibly.
#' @method print summary.act_BayesGLM_cifti
print.summary.act_BayesGLM_cifti <- function(x, ...) {
  cat("====BayesGLM_cifti Activations==============\n")
  cat(paste0("Activated locations:\n"))
  if (!("active" %in% names(x$activations[[1]]))) {
    # for (ii in seq(length(x$activations))) {
    #   names(x$activations[[ii]]) <- paste0(
    #     names(x$activations[[ii]]), ", ", names(x$activations)[ii]
    #   )
    # }
    x$activations <- do.call(c, x$activations)
  }
  for (ii in seq(length(x$activations))) {
    act_ii <- x$activations[[ii]]
    cat(paste0(
      "  ", names(x$activations)[ii],
      " (", sum(act_ii[,1]), " modeled locations)\n"
    ))
    for (kk in seq(ncol(act_ii))) {
      cat(paste0(
        "    ", colnames(act_ii)[kk],
        ": ", act_ii["TRUE",kk], "\n"
      ))
    }
  }
  cat("GLM type:   ", x$method, "\n")
  cat("alpha:      ", x$alpha, "\n")
  cat("gamma:      ", x$gamma, "\n")
  if (x$correction != "not applicable") {
    cat("Correction: ", x$correction, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.act_BayesGLM_cifti
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print act_BayesGLM_cifti
print.act_BayesGLM_cifti <- function(x, ...) {
  print.summary.act_BayesGLM_cifti(summary(x))
}
