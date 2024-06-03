#' Summarize a \code{"act_BGLM0"} object
#'
#' Summary method for class \code{"act_BGLM0"}
#'
#' @param object Object of class \code{"act_BGLM0"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.act_BGLM0"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary act_BGLM0
summary.act_BGLM0 <- function(object, ...) {
  # This is the possible structure of `object`:
  # object$activations$cortex_left$(session1)$`gamma=g`$active # CIFTI BayesGLM
  # object$activations$(session1)$`gamma=g`$active # BayesGLM

  # Try to make one consistent structure...
  act <- object$activations[!vapply(object$activations, is.null, FALSE)]
  if (inherits(object, "act_BGLM0")) { act <- list(model=act) } # See line near end of function
  if ("p_values" %in% names(act[[1]][[1]])) { act <- list(single_session=act) }
  act <- act[!vapply(act, function(q){is.null(q[[1]])}, FALSE)]

  # Now we have a common form:
  # act$cortex_left$session1$`gamma=g`$active # CIFTI BayesGLM
  # act$model$session1$`gamma=g`$active # BayesGLM

  x <- list(
    activations_g1 = lapply(act, function(x_model){
      lapply(x_model, function(x_sess){
        x_sess <- x_sess[[1]] # use first level of gamma
        q <- apply(x_sess$active, 2, function(avec){c(`TRUE`=sum(avec>0), `FALSE`=sum(avec==0))})
        colnames(q) <- object$field_names
        q
      })
    }),
    method=object$method,
    alpha=object$alpha,
    gamma=object$gamma,
    correction=object$correction
    #excur_method
  )

  if (inherits(object, "act_BGLM0")) { x$activations_g1 <- x$activations_g1$model } # See line near top of function

  class(x) <- "summary.act_BGLM0"
  x
}

#' @rdname summary.act_BGLM0
#' @export
#'
#' @param x Object of class \code{"summary.act_BGLM0"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.act_BGLM0
print.summary.act_BGLM0 <- function(x, ...) {
  gamma_msg <- if (length(x$gamma)>1) {
    paste0(" (lowest threshold, gamma=", x$gamma[1], ")")
  } else {
    ""
  }
  cat("====BayesGLM0 Activations===================\n")
  cat(paste0(
    "Activated locations (",
    sum(x$activations_g1[[1]][,1]), " modeled locations)",
    gamma_msg, ":\n"
  ))
  for (ii in seq(length(x$activations_g1))) {
    act_ii <- x$activations_g1[[ii]]
    cat(paste0(
      "  ", names(x$activations_g1)[ii],
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
  cat("gamma:      ", paste0(x$gamma, collapse=", "), "\n")
  if (x$correction != "not applicable") {
    cat("Correction: ", x$correction, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.act_BGLM0
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print act_BGLM0
print.act_BGLM0 <- function(x, ...) {
  print.summary.act_BGLM0(summary(x))
}

#' Summarize a \code{"act_BGLM"} object
#'
#' Summary method for class \code{"act_BGLM"}
#'
#' @param object Object of class \code{"act_BGLM"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A \code{"summary.act_BGLM"} object, a list summarizing the
#'  properties of \code{object}.
#' @method summary act_BGLM
summary.act_BGLM <- function(object, ...) {
  x <- summary.act_BGLM(object, ...)
  class(x) <- "summary.act_BGLM"
  x
}

#' @rdname summary.act_BGLM
#' @export
#'
#' @param x Object of class \code{"summary.act_BGLM"}.
#' @return \code{NULL}, invisibly.
#' @method print summary.act_BGLM
print.summary.act_BGLM <- function(x, ...) {
  gamma_msg <- if (length(x$gamma)>1) {
    paste0(" (lowest threshold, gamma=", x$gamma[1], ")")
  } else {
    ""
  }
  cat("====BayesGLM Activations====================\n")
  cat(paste0("Activated locations", gamma_msg, ":\n"))
  if (!("active" %in% names(x$activations_g1[[1]]))) {
    # for (ii in seq(length(x$activations_g1))) {
    #   names(x$activations_g1[[ii]]) <- paste0(
    #     names(x$activations_g1[[ii]]), ", ", names(x$activations_g1)[ii]
    #   )
    # }
    x$activations_g1 <- do.call(c, x$activations_g1)
  }
  for (ii in seq(length(x$activations_g1))) {
    act_ii <- x$activations_g1[[ii]]
    cat(paste0(
      "  ", names(x$activations_g1)[ii],
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
  cat("gamma:      ", paste0(x$gamma, collapse=", "), "\n")
  if (x$correction != "not applicable") {
    cat("Correction: ", x$correction, "\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.act_BGLM
#' @export
#'
#' @return \code{NULL}, invisibly.
#' @method print act_BGLM
print.act_BGLM <- function(x, ...) {
  print.summary.act_BGLM(summary(x))
}
