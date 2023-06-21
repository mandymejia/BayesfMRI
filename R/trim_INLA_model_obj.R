#' Trim INLA object
#' 
#' Trim an INLA object to only include what is necessary for 
#'  \code{id_activations} or \code{BayesGLM2}.
#'
#' @param INLA_model_obj An object of class \code{"inla"}.
#' @param minimal Just keep the two parameters needed for \code{BayesGLM2}?
#'  Default: \code{FALSE}. \code{!minimal} is required for 
#'  \code{id_activations}, but \code{minimal} is sufficient for 
#'  \code{BayesGLM2}.
#'
#' @return A trimmed \code{"inla"} object.
#' @keywords internal
trim_INLA_model_obj <- function(INLA_model_obj, minimal=FALSE) {
  if (!inherits(INLA_model_obj, "inla")) {
    stop("This function only applies to objects with the 'inla' class.")
  }

  out_object <- list()
  out_object$misc$theta.mode <- INLA_model_obj$misc$theta.mode
  out_object$misc$cov.intern <- INLA_model_obj$misc$cov.intern
  if (!minimal) {
    out_object$.args$control.compute$config <- INLA_model_obj$.args$control.compute$config
    out_object$marginals.random <- INLA_model_obj$marginals.random
    out_object$misc$configs <- INLA_model_obj$misc$configs
  }

  class(out_object) <- "inla"
  out_object
}
