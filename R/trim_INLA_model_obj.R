#' Trim an INLA object to only include what is necessary for the id_activations function
#'
#' @param INLA_model_obj an object of class 'inla' that will be trimmed
#'
#' @return a trimmed object of class 'inla'
#' @keywords internal
trim_INLA_model_obj <- function(INLA_model_obj) {
  if(!inherits(INLA_model_obj, "inla"))
    stop("This function only applies to objects with the 'inla' class.")
  out_object <- list(
    .args = list(
      control.compute = list(
        config = INLA_model_obj$.args$control.compute$config
      )
    ),
    marginals.random = INLA_model_obj$marginals.random,
    misc = list(
      configs = INLA_model_obj$misc$configs
    )
  )
  class(out_object) <- "inla"
  return(out_object)
}
