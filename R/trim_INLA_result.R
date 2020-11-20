#' Trim an INLA object to only include what is necessary for the id_activations function
#'
#' @param INLA_result an object of class 'inla' that will be trimmed
#'
#' @return a trimmed object of class 'inla'
#' @keywords internal
trim_INLA_result <- function(INLA_result) {
  if(class(INLA_result) != "inla")
    stop("This function only applies to objects with the 'inla' class.")
  out_object <- list(
    .args = list(
      control.compute = list(
        config = INLA_result$.args$control.compute$config
      )
    ),
    marginals.random = INLA_result$marginals.random,
    misc = list(
      configs = list(
        nconfig = INLA_result$misc$configs$nconfig,
        contents = INLA_result$misc$configs$contents,
        config = sapply(INLA_result$misc$configs$config, function(conf) {
          return(conf[c('mean','Q','Qinv','log.posterior','pred_idx')])
        }, simplify = F)
      )
    )
  )
  class(out_object) <- "inla"
  return(out_object)
}
