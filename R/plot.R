#' Plot a 2D image using ggplot
#'
#' @param X A matrix or list of matrices to be plotted
#' @param color_palette A color palette as a data frame with two columns. The
#'   first column should contain colors (in HEX, character, or other format),
#'   and the second column should contain locations for color breaks on the unit
#'   interval describing the locations for color switching. This will default to
#'   the output from \code{ROY_BIG_BL} palette from the \code{ciftiTools}
#'   package.
#' @param zlim A vector of length 2 describing the endpoints of the color
#'   palette.
#'
#' @return A ggplot graphical object.
#' 
#' @importFrom ciftiTools ROY_BIG_BL
#' @importFrom reshape2 melt
#' @importFrom utils data
#' 
#' @export
plot_slice <- function(X, color_palette = NULL, zlim = NULL) {

  # Check to see that the INLA package is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_slice` requires the `ggplot` package. Please install it.", call. = FALSE)
  }

  # Check to see that the INLA package is installed
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("`plot_slice` requires the `purrr` package. Please install it.", call. = FALSE)
  }

  if(class(X) == "matrix") X = list(single_activation_field = X)
  if(class(X) != "list") stop("Expected a matrix or list for X.")
  if(any(!sapply(X,function(x) {
    "matrix" %in% class(x)
  }, simplify = T))) {
    stop("All list images should be matrices.")
  }

  if(is.null(zlim)) {
    zmin <- min(reshape2::melt(X)$value,na.rm = T)
    zmax <- max(reshape2::melt(X)$value,na.rm = T)
    zlim <- c(zmin,zmax)
  }
  if(is.null(color_palette)){
    color_palette <- ciftiTools::ROY_BIG_BL(
      min = zlim[1],
      max = zlim[2],
      mid = mean(zlim),
      pos_half = FALSE
    )
  }
  if(min(reshape2::melt(X)$value, na.rm = T) >= 0) {
    color_palette <- ciftiTools::ROY_BIG_BL(
      min = zlim[1],
      max = zlim[2],
      mid = mean(zlim),
      pos_half = TRUE
    )
  }

  # Hacky way to avoid R CMD CHECK problems. The other solution is Importing ggplot2.
  ggplot <- geom_raster <- aes <- scale_fill_gradientn <- facet_grid <- NULL
  labs <- theme_bw <- theme <- element_blank <- NULL
  Var1 <- Var2 <- value <- NULL

  out_grob <- reshape2::melt(X)
  out_grob$value <- ifelse(out_grob$value < min(zlim, na.rm = T), min(zlim, na.rm = T), value)
  out_grob$value <- ifelse(out_grob$value > max(zlim, na.rm = T), max(zlim, na.rm = T), value)
  
  ggplot(out_grob) +
    geom_raster(aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_gradientn("",colors = rev(color_palette$color),
                         limits = zlim,
                         na.value = "white") +
    facet_grid(.~L1) +
    labs(x="", y="") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
}

#' Plot results from a BayesGLM object for 2D Analyses
#'
#' @param BayesGLM_object An object of class "BayesGLM"
#' @param session_name The name of the session to plot the results from (defaults to the first session)
#' @param zlim The color limits for plotting the coefficient values. Defaults to the minimum and maximum of the point estimates
#'
#' @return A ggplot2 object
#' 
#' @importFrom ciftiTools ROY_BIG_BL
#' @importFrom INLA inla.spde.make.A
#' 
#' @export
plot_BayesGLM_slice <- function(BayesGLM_object, session_name = NULL, zlim = NULL) {
  Var1 <- Var2 <- value <- NULL # to prevent package build warning at ggplot line

  # Check to see that the INLA package is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_slice` requires the `ggplot` package. Please install it.", call. = FALSE)
  }

  # Check to see that the INLA package is installed
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("`plot_slice` requires the `purrr` package. Please install it.", call. = FALSE)
  }


  # Create a conversion matrix
  in_binary_mask <- which(BayesGLM_object$mask == 1, arr.ind = T)
  in_binary_mask <- in_binary_mask[,2:1]
  convert_mat_A <- INLA::inla.spde.make.A(mesh = BayesGLM_object$mesh, loc = in_binary_mask)
  # Extract the point estimates
  if(is.null(session_name)) session_name <- BayesGLM_object$session_names
  point_estimates <- sapply(session_name, function(sn){
    as.matrix(convert_mat_A %*% BayesGLM_object$beta_estimates[[sn]])
  }, simplify = F)
  if(is.null(zlim)) zlim <- c(min(unlist(point_estimates)),
                              max(unlist(point_estimates)))
  wb_palette <- ciftiTools::ROY_BIG_BL(min = zlim[1], max = zlim[2], mid = mean(zlim), pos_half = FALSE)
  coef_images <- sapply(point_estimates, function(pe) {
    out <- sapply(split(pe, col(pe)), function(beta) {
      beta_out <- BayesGLM_object$mask
      beta_out[beta_out == 1] <- beta
      beta_out[beta_out == 0] <- NA
      return(beta_out)
    }, simplify = F)
    names(out) <- BayesGLM_object$beta_names
    return(out)
  }, simplify= F)

  reshape2::melt(coef_images) %>%
    ggplot() +
    geom_raster(aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_gradientn("",colors = rev(wb_palette$color),
                         # values = wb_palette$value,
                         limits = zlim,
                         na.value = "white") +
    facet_grid(L1~L2) +
    labs(x="", y="") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
}

#' S3 method: use \code{plot_BayesGLM_slice} to plot a \code{"BayesGLM"} object
#' 
#' @param x An object of class "BayesGLM"
#' @param ... Additional arguments to \code{\link{plot_BayesGLM_slice}}
#'
#' @method plot BayesGLM
#' 
#' @export
#' 
plot.BayesGLM <- function(x, ...){
  plot_BayesGLM_slice(x, ...)
}