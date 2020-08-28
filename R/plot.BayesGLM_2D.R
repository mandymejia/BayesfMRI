#' Plot results from a BayesfMRI object for 2D Analyses
#'
#' @param BayesGLM_object An object of class "BayesfMRI"
#' @param mask The mask used to convert from cortical surface to volumetric data
#' @param session_name The name of the session to plot the results from (defaults to the first session)
#' @param zlim The color limits for plotting the coefficient values. Defaults to the minimum and maximum of the point estimates
#'
#' @return A ggplot2 object
#' @importFrom ciftiTools ROY_BIG_BL
#' @import dplyr
#' @import ggplot2
#' @export
plot.BayesGLM_2d <- function(BayesGLM_object, mask, session_name = NULL, zlim = NULL) {
  # Create a conversion matrix
  in_binary_mask <- which(mask == 1, arr.ind = T)
  in_binary_mask <- in_binary_mask[,2:1]
  convert_mat_A <- inla.spde.make.A(mesh = BayesGLM_object$mesh, loc = in_binary_mask)
  # Extract the point estimates
  if(is.null(session_name)) session_name <- BayesGLM_object$session_names
  point_estimates <- sapply(session_name, function(sn){
    as.matrix(convert_mat_A %*% BayesGLM_object$beta_estimates[[sn]])
  }, simplify = F)
  if(is.null(zlim)) zlim <- c(min(unlist(point_estimates)),
                              max(unlist(point_estimates)))
  div_pal <- ROY_BIG_BL(min = zlim[1], max = zlim[2],mid = 0, pos_half = FALSE)
  coef_images <- sapply(point_estimates, function(pe) {
    out <- sapply(split(pe, col(pe)), function(beta) {
      beta_out <- mask
      beta_out[beta_out == 1] <- beta
      beta_out[beta_out == 0] <- NA
      return(beta_out)
    }, simplify = F)
    names(out) <- BayesGLM_object$beta_names
    return(out)
  }, simplify= F)
  out_grob <- reshape2::melt(coef_images) %>%
    ggplot() +
    geom_raster(aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_gradient2("") +
    facet_grid(L1~L2) +
    labs(x="", y="") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  return(out_grob)
}
