#' Plot results from a BayesfMRI object
#'
#' @param BayesGLM_object An object of class "BayesfMRI"
#' @param mask The mask used to convert from cortical surface to volumetric data
#' @param session_name The name of the session to plot the results from (defaults to the first session)
#'
#' @return Currently a list of final coefficient array images. This will change in future versions.
#' @keywords internal
plot.BayesGLM <- function(BayesGLM_object, mask, session_name = NULL) {
  # Create a conversion matrix
  in_binary_mask <- which(mask == 1, arr.ind = T)
  in_binary_mask <- in_binary_mask[,2:1]
  convert_mat_A <- inla.spde.make.A(mesh = BayesGLM_object$mesh, loc = in_binary_mask)
  # Extract the point estimates
  if(is.null(session_name)) session_name <- BayesGLM_object$session_names[1]
  point_estimates <- as.matrix(convert_mat_A %*% BayesGLM_object$beta_estimates[[session_name]])
  # Create the image
  coef_images <- sapply(split(point_estimates, col(point_estimates)), function(beta) {
    beta_out <- mask
    beta_out[beta_out == 1] <- beta
    beta_out[beta_out == 0] <- NA
    return(beta_out)
  }, simplify = F)
  reshape2::melt(coef_images) %>%
    ggplot() +
    geom_raster(aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_gradient2("") +
    facet_grid(.~L1)
  # Using Mandy's code
  plot(mesh,rgl=TRUE, draw.edges=FALSE)
  nColors <- 64
  colindex <- as.integer(cut(BayesGLM_object$beta_estimates$single_session[,1],breaks=nColors))
  plot(mesh,col=heat.colors(nColors)[colindex], rgl=TRUE)
  return(coef_images)
}
