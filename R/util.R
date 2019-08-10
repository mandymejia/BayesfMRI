
#' Summarise BayesGLMfMRI objects
#'
#' Summary method for class "BayesGLMfMRIobj"
#'
#' @param object an object of class "BayesGLMfMRIobj"
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary BayesGLMfMRIobj
summary.BayesGLMfMRIobj <- function(object, ...)
{
  out <- list()
  class(out) <- "summary.BayesGLMfMRIobj"
  out$sessions = object$sessions
  out$call <- object$call
  out$inla.summary <- summary(object$model)
  return(out)
}


#' @param x an object of class "summary.BayesGLMfMRIobj"
#' @export
#' @method print summary.BayesGLMfMRIobj
#' @rdname summary.BayesGLMfMRIobj
print.summary.BayesGLMfMRIobj <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("Sessions: ", x$sessions,"\n")
  cat("Time used:\n", x$inla.summary$cpu.used)
}

#' @export
#' @method print BayesGLMfMRIobj
#' @rdname summary.BayesGLMfMRIobj
print.BayesGLMfMRIobj <- function(x, ...) {
  print.summary.BayesGLMfMRIobj(summary(x))
}
