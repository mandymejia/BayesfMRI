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

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_slice` requires the `ggplot` package. Please install it.", call. = FALSE)
  }

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
  if(min(unlist(X), na.rm = T) >= 0) {
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

  X_df <- reshape2::melt(X)
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
#' @importFrom INLA inla.spde.make.A
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


  # ggplot(reshape2::melt(coef_images)) +
  point_df <- reshape2::melt(point_estimates)
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

#' S3 method: use \code{\link[ciftiTools]{view_xifti_surface}} to plot a \code{"BayesGLM_CIFTI"} object
#'
#' @param x An object of class "BayesGLM_CIFTI"
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
#' @method plot BayesGLM_CIFTI
#'
# @importFrom ciftiTools view_xifti_surface
#' @export
#'
plot.BayesGLM_CIFTI <- function(x, session=NULL, method=NULL, idx=NULL, zlim=c(-1, 1), ...){

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
