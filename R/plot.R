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
#' @importFrom utils data
#'
#' @export
plot_slice <- function(X, color_palette = NULL, zlim = NULL) {
  # Hacky way to avoid R CMD CHECK problems. The other solution is Importing ggplot2.
  ggplot <- geom_raster <- aes <- scale_fill_gradientn <- facet_grid <- NULL
  labs <- theme_bw <- theme <- element_blank <- NULL
  Var1 <- Var2 <- value <- NULL
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_slice` requires the `ggplot` package. Please install it.", call. = FALSE)
  }

  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("`plot_slice` requires the `purrr` package. Please install it.", call. = FALSE)
  }

  if(inherits(X, "matrix")) X = list(single_activation_field = X)
  if(!inherits(X, "list")) stop("Expected a matrix or list for X.")
  if (!all(sapply(X, inherits, "matrix", simplify=TRUE))) {
    stop("All list images should be matrices.")
  }

  if(is.null(zlim)) {
    zmin <- min(melt_mat2(X)$value,na.rm = T)
    zmax <- max(melt_mat2(X)$value,na.rm = T)
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
  if(min(unlist(X), na.rm = T) >= 0) {
    color_palette <- ciftiTools::ROY_BIG_BL(
      min = zlim[1],
      max = zlim[2],
      mid = mean(zlim),
      pos_half = TRUE
    )
  }



  X_df <- melt_mat2(X)
  X_df$value <- ifelse(X_df$value < min(zlim, na.rm = T), min(zlim, na.rm = T), X_df$value)
  X_df$value <- ifelse(X_df$value > max(zlim, na.rm = T), max(zlim, na.rm = T), X_df$value)

  out_grob <- ggplot(X_df) +
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

#' S3 method: use \code{\link[ciftiTools]{view_xifti_surface}} to plot a \code{"BayesGLM_cifti"} object
#'
#' @param x An object of class "BayesGLM_cifti"
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param method "Bayesian" or "classical". \code{NULL} (default) will use
#'  the Bayesian results if available, and the classical results if not.
#' @param idx The data columns to plot. Overrides the \code{idx} argument to
#'  \code{\link[ciftiTools]{view_xifti_surface}}. \code{NULL} (default) will
#'  use all the data columns.
#' @param zlim Overrides the \code{zlim} argument for
#'  \code{\link[ciftiTools]{view_xifti_surface}}. Default: \code{c(-1, 1)}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti_surface}}
#'
#' @method plot BayesGLM_cifti
#'
# @importFrom ciftiTools view_xifti_surface
#' @export
#'
plot.BayesGLM_cifti <- function(x, session=NULL, method=NULL, idx=NULL, zlim=c(-1, 1), ...){

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("This function requires the `ciftiTools` package. Please install it.")
  }

  # Method
  if (is.null(method)) {
    method <- ifelse(
      is.null(x$betas_Bayesian) || is.null(x$betas_Bayesian[[1]]), 
      "classical", "Bayesian"
    )
  }
  method <- match.arg(method, c("classical", "Bayesian"))
  method <- paste0("betas_", method)
  if (is.null(x[[method]])) {
    stop(paste("Method", gsub("betas_", "", method, fixed=TRUE), "does not exist."))
  }

  # Session
  if (is.null(session)) { 
    if (length(x[[method]]) > 1) { message("Plotting the first session.") }
    session <- 1 
  }
  if (is.null(x[[method]][[session]])) {
    stop(paste("Session", session, "of method", method, "does not exist."))
  }

  # Column index
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, x[[method]][[session]]$data)))
  }

  # Plot
  ciftiTools::view_xifti_surface(x[[method]][[session]], idx=idx, zlim=zlim, ...)
}
