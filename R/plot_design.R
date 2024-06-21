#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"BGLM"} object
#'
#' @param x An object of class "BfMRI_design".
#' @param ... Additional arguments to \code{\link{plot_design}}.
#' @method plot BfMRI_design
#' @export
#'
#' @return Result of the call to \code{\link{plot_design}}
#'
plot.BfMRI_design <- function(x, ...){
  stopifnot(inherits(x, "BfMRI_design"))
  plot_design(x$design, ...)
}

#' Plot design matrix
#'
#' @param design The timepoints by fields design matrix or data.frame.
#' @param method \code{"lineplot"} (default) or \code{"imageplot"}.
#' @param ... Additional arguments to \code{plot_design_line} or
#' \code{plot_design_image}.
#' @return A ggplot
#' @export
#'
plot_design <- function(design, method=c("lineplot", "imageplot"), ...){
  design <- as.data.frame(design)
  method <- match.arg(method, c("lineplot", "imageplot"))

  switch(method,
    lineplot=plot_design_line(design, ...),
    imageplot=plot_design_image(design, ...)
  )
}

#' Plot design with lineplot
#'
#' @rdname plot_design
#' @param design The timepoints by fields design matrix or data.frame.
#' @param colors The name of a ColorBrewer palette (see
#'  RColorBrewer::brewer.pal.info and colorbrewer2.org), the name of a
#'  viridisLite palette, or a character vector of colors. Default:
#'  \code{"Set1"}.
#' @param linetype,linewidth,alpha Parameters for \code{ggplot2::geom_line}. 
#'  Defaults: \code{"solid"} linetype, \code{0.7} linewidth and \code{0.8} 
#'  alpha. \code{linetype} can also be a vector of options with length matching
#'  the number of fields in \code{design}.
#' @return A ggplot
#' @export
#' @importFrom ciftiTools make_color_pal
#'
plot_design_line <- function(design, colors="Set1", linetype="solid", linewidth=.7, alpha=.8){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please download the `ggplot2` package.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Please download the `tidyr` package.")
  }

  df <- as.data.frame(design)
  nT <- nrow(df)
  nK <- ncol(df)

  colors <- suppressWarnings(
    make_color_pal(colors, "qualitative", zlim=nK)$color
  )

  df <- cbind(df, data.frame(idx=seq(nT)))
  df <- tidyr::pivot_longer(df, seq(nK))
  colnames(df)[colnames(df)=="name"] <- "Field"
  df$Field <- factor(df$Field, levels=colnames(design))

  if (length(linetype) == nK) {
    plt <- ggplot2::ggplot(df, ggplot2::aes_string(x="idx", y="value", col="Field", linetype="Field")) +
      ggplot2::geom_hline(yintercept=0, color="black", linetype="dashed") +
      ggplot2::geom_line(linewidth=linewidth, alpha=alpha) +
      ggplot2::scale_linetype_manual(values=linetype)
  } else {
    plt <- ggplot2::ggplot(df, ggplot2::aes_string(x="idx", y="value", col="Field")) +
      ggplot2::geom_hline(yintercept=0, color="black", linetype="dashed") +
      ggplot2::geom_line(linetype=linetype, linewidth=linewidth, alpha=alpha)
  }

  plt + ggplot2::theme_bw() +
    ggplot2::xlab("Volume") + ggplot2::ylab("Value") +
    ggplot2::scale_color_manual(values=colors) +
    ggplot2::scale_x_continuous(limits=c(0-1e-8, nT+1e-8), expand=c(0,0))
}

#' Plot design with imageplot
#'
#' @rdname plot_design
#' @param design The timepoints by fields design matrix or data.frame.
#' @return A ggplot
#' @export
#'
plot_design_image <- function(design){
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Please download the `grDevices` package.")
  }
  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop("Please download the `graphics` package.")
  }

  design <- as.matrix(design)
  nT <- nrow(design)
  nK <- ncol(design)

  graphics::image(
    t(design[seq(nT,1),]), # rev
    xlab="Field", ylab="Volume", col=grDevices::gray.colors(256), axes=FALSE
  )
  k_divs <- if (nK == 1) { 0 } else { seq(0, nK-1)/(nK-1) }
  graphics::axis(1, at = k_divs, labels=colnames(design))
  t_divs <-seq(0, 3)/3
  graphics::axis(2, at = t_divs, labels=rev(round(t_divs*nT))) # rev
}
