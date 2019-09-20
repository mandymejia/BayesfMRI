#' Extracts posterior density estimates for hyperparameters
#'
#' @param object An object of class ‘"inla"’, a result of a call to \code{inla()}
#' @param spde The model used for the latent fields in the \code{inla()} call, an object of class ‘"inla.spde"’
#'
#' @return Long-form data frame containing posterior densities for the hyperparameters associated with each latent field
#' @export
#' @importFrom INLA inla.spde2.result
#'
#' @examples \dontrun{}
get_posterior_densities_vol3D <- function(object, spde){

  beta_names <- names(object$summary.random)

  for(b in beta_names){
    df.b <- inla.extract.el(object$marginals.hyperpar,
                                     paste("Theta[^ ]+ for ", b, "$", sep = ""))
    df.b <- as.data.frame(df.b)
    #colnames(df.b) <-

    df.b$beta <- b
    if(b==beta_names[1]) df <- df.b else df <- cbind(df, df.b)
  }
  # df <- df[,c('beta','param','value','density')]
  return(df)
}

