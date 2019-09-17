
#' Summarise BayesGLM objects
#'
#' Summary method for class "BayesGLMobj"
#'
#' @param object an object of class "BayesGLMobj"
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary BayesGLMobj
summary.BayesGLMobj <- function(object, ...)
{
  out <- list()
  class(out) <- "summary.BayesGLMobj"
  out$sessions = object$sessions
  out$call <- object$call
  out$inla.summary <- summary(object$model)
  return(out)
}


#' @param x an object of class "summary.BayesGLMobj"
#' @export
#' @method print summary.BayesGLMobj
#' @rdname summary.BayesGLMobj
print.summary.BayesGLMobj <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("Sessions: ", x$sessions,"\n")
  cat("Time used:\n", x$inla.summary$cpu.used)
}

#' @export
#' @method print BayesGLMobj
#' @rdname summary.BayesGLMobj
print.BayesGLMobj <- function(x, ...) {
  print.summary.BayesGLMobj(summary(x))
}

# TO DO: Add print and summary functions for session object (may need as.session function too, and move is.session here)
