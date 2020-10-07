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
#' @import ggplot2
#' @import purrr
#' @importFrom ciftiTools ROY_BIG_BL
#' @importFrom reshape2 melt
#' @importFrom utils data
#' @importFrom dplyr mutate
#' @export
plot_slice <- function(X, color_palette = NULL, zlim = NULL) {
  if(class(X) == "matrix") X = list(single_activation_field = X)
  if(class(X) != "list") stop("Expected a matrix or list for X.")
  if(any(!sapply(X,function(x) {
    "matrix" %in% class(x)
  }, simplify = T))) {
    stop("All list images should be matrices.")
  }
  # Make R CMD check happy
  Var1 <- Var2 <- value <- NULL
  requireNamespace("ggplot2")
  requireNamespace("purrr")
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
  out_grob <- reshape2::melt(X) %>%
    dplyr::mutate(value = ifelse(value < min(zlim,na.rm = T),min(zlim,na.rm = T),value),
                  value = ifelse(value > max(zlim,na.rm = T),max(zlim,na.rm = T),value)) %>%
    ggplot() +
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
