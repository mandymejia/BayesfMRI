#' Extracts posterior density estimates for hyperparameters
#'
#' @param object An object of class ‘"inla"’, a result of a call to \code{inla()}
#' @param spde The model used for the latent fields in the \code{inla()} call, an object of class ‘"inla.spde"’
#'
#' @return Long-form data frame containing posterior densities for the hyperparameters associated with each latent field
#' @export
#' @importFrom INLA inla.spde2.result
#' @importFrom INLA inla.extract.el
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
get_posterior_densities_vol3D <- function(object, spde){

  hyper_names <- names(object$marginals.hyperpar)

  for(h in hyper_names){
    df.h <- inla.extract.el(object$marginals.hyperpar, h)
    df.h <- as.data.frame(df.h)
    names(df.h) <- c('x','y')
    df.h$hyper_name <- h
    df.h$beta_name <- ifelse(grepl('bbeta',h), #if bbeta appears in name
                        gsub('.+bbeta','bbeta',h), #rename as bbeta*
                        NA) #NA if this is the precision hyperparameter
    df.h$theta_name <- ifelse(grepl('Theta',h), #if bbeta appears in name
                              gsub(' for.+','',h), #rename as Theta*
                              NA) #NA if this is the precision hyperparameter

    if(h==hyper_names[1]) df <- df.h else df <- rbind(df, df.h)
  }
  return(df)
}

