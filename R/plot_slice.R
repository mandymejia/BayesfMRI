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

  # Hacky way to avoid R CMD CHECK problems. The other solution is Importing ggplot2.
  ggplot <- geom_raster <- aes <- scale_fill_gradientn <- facet_grid <- NULL
  labs <- theme_bw <- theme <- element_blank <- NULL

  if(class(X) == "matrix") X = list(single_activation_field = X)
  if(class(X) != "list") stop("Expected a matrix or list for X.")
  if(any(!sapply(X,function(x) {
    "matrix" %in% class(x)
  }, simplify = T))) {
    stop("All list images should be matrices.")
  }
  # Make R CMD check happy
  Var1 <- Var2 <- value <- NULL
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
  out_grob <- reshape2::melt(X)
  out_grob$value <- ifelse(out_grob$value < min(zlim, na.rm = T), min(zlim, na.rm = T), value)
  out_grob$value <- ifelse(out_grob$value > max(zlim, na.rm = T), max(zlim, na.rm = T), value)
  out_grob <- ggplot(out_grob) +
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
  return(out_grob)
}
