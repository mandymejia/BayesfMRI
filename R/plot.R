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

#' Plot results from a BayesGLM object for 2D Analyses
#'
#' @param BayesGLM_object An object of class "BayesGLM"
#' @param session_name The name of the session to plot the results from (defaults to the first session)
#' @param zlim The color limits for plotting the coefficient values. Defaults to the minimum and maximum of the point estimates
#'
#' @return A ggplot2 object
#'
#' @importFrom ciftiTools ROY_BIG_BL
#'
#' @export
plot_BayesGLM_slice <- function(BayesGLM_object, session_name = NULL, zlim = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_slice` requires the `ggplot2` package. Please install it.", call. = FALSE)
  }
  # to prevent package build warning at ggplot lines
  ggplot <- geom_raster <- aes <- scale_fill_gradientn <- facet_grid <- NULL
  labs <- theme_bw <- theme <- element_blank <- Var1 <- Var2 <- value <- NULL

  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("`plot_slice` requires the `purrr` package. Please install it.", call. = FALSE)
  }

  check_INLA(FALSE)

  # Create a conversion matrix
  in_binary_mask <- which(BayesGLM_object$mask == 1, arr.ind = T)
  in_binary_mask <- in_binary_mask[,2:1]
  mesh <- make_slice_mesh(BayesGLM_object$mask)
  convert_mat_A <- INLA::inla.spde.make.A(mesh = mesh, loc = in_binary_mask)
  # Extract the point estimates
  if(is.null(session_name)) session_name <- BayesGLM_object$GLMs_Bayesian$session_names
  # point_estimates <- sapply(session_name, function(sn){
  #   as.matrix(convert_mat_A %*% BayesGLM_object$beta_estimates[[sn]])
  # }, simplify = F)
  point_estimates <- BayesGLM_object$betas_Bayesian
  if(is.null(zlim)) zlim <- c(min(unlist(point_estimates), na.rm = T),
                              max(unlist(point_estimates), na.rm = T))
  wb_palette <- ciftiTools::ROY_BIG_BL(min = zlim[1], max = zlim[2], mid = mean(zlim), pos_half = FALSE)
  # coef_images <- sapply(point_estimates, function(pe) {
  #   out <- sapply(split(pe, col(pe)), function(beta) {
  #     beta_out <- BayesGLM_object$mask
  #     beta_out[beta_out == 1] <- beta
  #     beta_out[beta_out == 0] <- NA
  #     return(beta_out)
  #   }, simplify = F)
  #   names(out) <- BayesGLM_object$beta_names
  #   return(out)
  # }, simplify= F)


  # ggplot(melt_mat2(coef_images)) +
  point_df <- melt_mat2(point_estimates)
  if(!is.null(BayesGLM_object$GLMs_Bayesian)) {
    point_df$L2 <- BayesGLM_object$GLMs_Bayesian$beta_names[point_df$L2]
  }
  out_grob <- ggplot(point_df) +
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
  return(out_grob)
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
    method <- ifelse(is.null(x$betas_Bayesian), "classical", "Bayesian")
  }
  method <- match.arg(method, c("classical", "Bayesian"))
  method <- paste0("betas_", method)
  if (is.null(x[[method]])) {
    stop(paste("Method", gsub("betas_", "", method, fixed=TRUE), "does not exist."))
  }

  # Session
  if (is.null(session)) { session <- 1 }
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

#' Melt a matrix into a data frame with identifiers for row and column
#'
#' This is meant to be a lightweight replacement for \code{reshape2::melt} that
#' will only work with matrices.
#'
#' @param x A matrix
#'
#' @return A data frame with the number of rows equal to the number of elements
#'   in \code{x}.
#' @export
#'
# @examples
# x <- matrix(rnorm(9),3,3)
# melt_mat(x)
melt_mat <- function(x) {
  if (!inherits(x, "matrix")) stop("x must have the matrix class.")
  out <- data.frame(row = c(row(x)), col = c(col(x)), value = c(x))
  return(out)
}

#' Melt list of matrices
#' 
#' See \code{melt_mat}
#' 
#' @param X_list list of matrices
#' @return A data.frame
#' @keywords internal
#' 
melt_mat2 <- function(X_list){
  X_list <- lapply(X_list, melt_mat)
  for (ii in seq(length(X_list))) {
    X_list[[ii]]$L1 = names(X_list)[ii]
  }
  out <- do.call(rbind, X_list)
  rownames(out) <- NULL
  out
}

#' Lightweight tile plot function
#'
#' This function has no dependencies outside of base R, but provides a
#'   reasonable approximation to the functionality of
#'   \code{ggplot2::geom_raster}.
#'
#' @param tile_df A data frame with three columns: \code{row}, \code{col}, and
#'   \code{value} describing locations within a matrix. See
#'   \code{\link{melt_mat}} for an example.
#' @param col A color palette
#' @param ncols The number of colors for a color palette, if \code{col} is not
#'   provided.
#' @param main Plot title (character)
#' @param na.color The color that should be used to represent \code{NA} values
#'   in the tile plot
#' @param zlim Color limits
#'
#' @return A tile plot done in base R graphics
#' @export
#' @importFrom grDevices heat.colors
#' @importFrom graphics axis layout par rect text
#'
# @examples
# x <- matrix(rnorm(50*50),50,50)
# x_df <- melt_mat(x)
# tile.plot(x_df)
tile.plot <- function(tile_df, col = NULL, ncols = NULL,
                      main = "", zlim = NULL, na.color  = "grey80") {
  if(inherits(tile_df, "matrix")) tile_df <- melt_mat(tile_df)
  .pardefault <- par()
  if(!is.null(col) & !is.null(ncols)) {
    warning("Defining ncols based on col.")
    ncols <- length(col)
  }
  if(is.null(col) & is.null(ncols)) {
    ncols <- 100
    col <- heat.colors(ncols)
  } else {
    if(is.null(ncols)) {
      ncols <- length(col)
    } else {
      col <- heat.colors(ncols)
    }
  }
  if(is.null(zlim)) zlim <- c(min(tile_df$value, na.rm = T), max(tile_df$value, na.rm = T))
  legend_ticks <- round(seq(zlim[1],zlim[2], length.out = 6),
                      digits = 3)
  prob_breaks <- seq(0,1,length.out = ncols)
  pb_diff <- prob_breaks[2] - prob_breaks[1]
  col_quants <- quantile(tile_df$value,na.rm = T,probs = prob_breaks)
  color_breaks <- quantile(seq(zlim[1],zlim[2],length.out = ncols),probs = prob_breaks)
  tile_cols <- vector("numeric",length(tile_df$value))
  tile_cols[tile_df$value >= color_breaks[ncols]] <- col[ncols]
  for(q in rev(seq(ncols))) {
    tile_cols[tile_df$value <= color_breaks[q]] <- col[q]
  }
  tile_cols[tile_cols == "0"] <- na.color
  rows <- max(tile_df$row)
  cols <- max(tile_df$col)
  cb_prime <- min(diff(color_breaks))
  par(mfrow = c(1,2), mar = c(1,1,2,1))
  layout(mat = matrix(c(1,2),nrow = 1, ncol = 2),
         widths = c(
           1.7, # tile plot width
           # max(c(1.7,dev.size()[1]*.85)), # tile plot width
           0.3 # legend width
           # min(c(0.3,dev.size()[1]*.15))
           )
         )
  plot(c(0,rows), c(0,cols), type = 'n', xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", main = main)
  rect(xleft = tile_df$row - 1,ybottom = tile_df$col - 1,
       xright = tile_df$row, ytop = tile_df$col, col = tile_cols,
       border = NA)
  par(mar=c(1,1,2,4))
  plot(c(0,1),c(zlim[1],zlim[2] + cb_prime), type = "n", xaxt = "n", yaxt = "n", xlab = "",
       ylab = "", bty = "n")
  rect(xleft = 0,
       ybottom = color_breaks,
       xright = 1,
       ytop = 1 + cb_prime,
       col = col, border = NA)
  axis(side = 4,at = c(legend_ticks[1],legend_ticks[-1]*(1+cb_prime)),
       labels = rep("",6), srt = 45, tck = 0.5)
  text(x = 1, adj = c(-1,0), pos = 4, y = c(legend_ticks[1],legend_ticks[-1]*(1+cb_prime)),
       labels = legend_ticks, srt = 0, xpd = NA)
  par(mfrow = .pardefault$mfrow, mar = .pardefault$mar)
  # suppressWarnings(par(.pardefault), classes = "warning")
}

