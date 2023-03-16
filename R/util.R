# #' Plot BayesfMRI.spde objects
# #'
# #' @param object Object of class BayesfMRI.spde (see \code{help(create_spde_vol3D)})
# #' @param colors (Optional) Vector of colors to represent each region.
# #' @param alpha Transparency level.
# #'
# #' @return
# #' @export
# #' @importFrom geometry tetramesh
# #' @import rgl
# #' @importFrom viridis viridis_pal
#
# plot.BayesfMRI.spde <- function(object, colors=NULL, alpha=0.5){
#   if(class(object) != 'BayesfMRI.spde') stop('object argument must be a BayesfMRI.spde object. See help(create_spde_vol3D).')
#   num_regions <- length(object$vertices)
#   if(is.null(colors)) colors <- viridis_pal()(num_regions)
#   if(length(colors) < num_regions) {
#     colors <- rep(colors, length.out=num_regions)
#     warning('Fewer colors specified than number of regions in the mesh.  Recycling colors to equal number of regions.')
#   }
#   for(ii in 1:num_regions){
#     if(ii==1) tetramesh(object$faces[[ii]], object$vertices[[ii]], col=colors[ii], alpha=alpha)
#     if(ii > 1) tetramesh(object$faces[[ii]], object$vertices[[ii]], clear=FALSE, col=colors[ii], alpha=alpha)
#
#   }
# }

# #' Plot BayesGLM objects
# #'
# #' Summary method for class "BayesGLM"
# #'
# #' @param object an object of class "BayesGLM"
# #' @param session_name NULL if BayesGLM object contains a single session; otherwise, the name of the session whose estimates to plot
# #' @param pal If NULL, viridis palette with 64 colors will be used.  Otherwise, specify a vector of color names.
# #' @param ... further arguments passed to or from other methods.
# #' @export
# #' @import viridis
# #' @method plot BayesGLM
# plot.BayesGLM <- function(object, session_name=NULL, pal=NULL, ...)
# {
#   session_names <- names(object$beta_estimates)
#
#   if((is.null(session_name)) & (length(session_names) > 1)) stop('If BayesGLM object includes multiple sessions, you must specify which session to plot.')
#   if(!is.null(session_name) & !(session_name %in% session_names)) stop('I expect the session_names argument to be one of the session names of the BayesGLM object, but it is not.')
#
#   if(is.null(session_name) & (length(session_names) == 1)) session_name <- session_names
#
#   ind <- which(names(object$beta_estimates) == session_name) #which element of list
#   est <- (object$beta_estimates)[[ind]]
#   K <- ncol(est)
#
#
#   if(is.null(pal)) {
#     nColors <- 64
#     pal <- viridis_pal(option='plasma', direction=1)(nColors)
#   } else {
#     if(min(areColors(pal)) < 1) stop('At least one of the elements of the pal argument is not a valid color representation.  See help(areColors).')
#     nColors <- length(pal)
#   }
#
#
#   for(k in 1:K){
#     x = est[,k]
#     colindex <- as.integer(cut(x,breaks=nColors))
#
#     #NEED TO CHECK WHICH TYPE OF BAYESGLM OBJECT (VOL OR CORTICAL) -- maybe use the mesh class?  or the spde_obj class?
#     #plot(mesh_LH_s$mesh,col=pal[colindex], rgl=TRUE)
#
#     #tetramesh(object$spde_obj$faces, object$spde_obj$vertices, col=pal[colindex], clear=FALSE)
#
#   }
#
#
# }
#

# #' Check whether each element of vector x is a valid color representation
# #'
# #' @param x Character vector
# #'
# #' @return A logical vector indicating which of the elements of x are valid color representations
# #' @importFrom grDevices col2rgb
# #' @export
#
# areColors <- function(x) {
#   sapply(x, function(X) {
#     tryCatch(is.matrix(col2rgb(X)),
#              error = function(e) FALSE)
#   })
# }


# #' Summarise BayesGLM objects
# #'
# #' Summary method for class "BayesGLM"
# #'
# #' @param object an object of class "BayesGLM"
# #' @param ... further arguments passed to or from other methods.
# #' @export
# #' @method summary BayesGLM
# summary.BayesGLM <- function(object, ...)
# {
#   out <- list()
#   class(out) <- "summary.BayesGLM"
#   out$sessions <- object$session_names
#   out$betas <- object$task_names
#   out$call <- object$INLA_result$call
#   out$inla.summary <- summary(object$model)
#   return(out)
# }


# #' @param x an object of class "summary.BayesGLM"
# #' @param ... further arguments passed to or from other methods.
# #' @export
# #' @method print summary.BayesGLM
# #' @rdname summary.BayesGLM
# print.summary.BayesGLM <- function(x, ...)
# {
#   cat("Call:\n")
#   print(x$call)
#   cat("Sessions: ", x$sessions,"\n")
#   cat("Time used:\n", x$inla.summary$cpu.used)
# }

# #' @export
# #' @method print BayesGLM
# #' @rdname summary.BayesGLM
# print.BayesGLM <- function(x, ...) {
#   print.summary.BayesGLM(summary(x))
# }

# TO DO: Add print and summary functions for session object (may need as.session function too, and move is.session here)



#' Find nonzero elements in matrix
#'
#' Find nonzero element in a matrix using 2-means clustering
#'
#' @importFrom stats kmeans
#'
#' @param beta_est A vector or matrix of values from which values close to zero
#'  should be assigned a value of zero.
#'
#' @return A vector or matrix of the same dimension as beta_est in which values
#'  close to zero are assigned the value of zero. The closeness of a value to
#'  zero is found by performing two-means clustering on the absolute values of
#'  beta_est, and ...
#'
#' @export
#'
find_nonzero <- function(beta_est) {
  vector_beta <- c(beta_est)
  if(any(is.na(vector_beta))) vector_beta <- vector_beta[!is.na(vector_beta)]
  km_beta <- kmeans(abs(vector_beta),2)
  which_nonzero <- which.max(km_beta$centers[,1])
  keep_nonzero <- as.numeric(km_beta$cluster == which_nonzero)
  out <- beta_est
  out[!is.na(out)] <- out[!is.na(out)] * keep_nonzero
  return(out)
}

#' Sequential 2-means variable selection
#'
#' @param x A vector consisting of all variables of interest for a single draw
#'  from a posterior distribution
#' @param b A scale parameter used to determine at what distance cluster centers
#'  are considered to be the same.
#'
#' @return The number of nonzero values detected within x
#'
#' @export
s2m <- function(x,b){
  two_means <- kmeans(abs(x),2)
  zero_idx <- which(two_means$cluster == which.min(two_means$centers))
  A <- x[zero_idx]
  two_centers <- kmeans(abs(A),2,algorithm=c("Lloyd"))
  iterations <- 1
  while(abs(two_centers$centers[1, 1] - two_centers$centers[2, 1]) > b) {
    zero_idx <- which(two_centers$cluster == which.min(two_centers$centers))
    A <- A[zero_idx]
    two_centers <- kmeans(abs(A),2,algorithm=c("Lloyd"))
    iterations <- iterations + 1
  }
  num_nonzero <- length(x) - length(A)
  return(num_nonzero)
}

#' Sequential 2-means on array B
#'
#' @param B An array of posterior samples (typically a matrix), in which the last margin corresponds to a single posterior sample
#' @param sigma  A scale parameter used to determine at what distance cluster centers are considered to be the same.
#'
#' @return An array of dimension `head(dim(B),-1)` with a point estimate of B based on the sequential 2-means method
#'
#' @importFrom stats quantile median
#'
#' @export
#'
#' @md
s2m_B <- function(B,sigma){
  nonzero_nums <- sapply(asplit(B,length(dim(B))),function(B_s) s2m(c(B_s),sigma))
  num_nonzero <- ceiling(median(nonzero_nums))
  median_B <- apply(B,seq(length(dim(B)) - 1),median)
  cutoff <- quantile(c(abs(median_B)),1 - num_nonzero/length(median_B))
  out <- median_B
  out[which(out < cutoff)] <- 0
  return(out)
}

#' Create a mask based on vertices that are invalid
#'
#' @param data A list of sessions, where each session is a list with elements
#'  BOLD, design and nuisance. See \code{?create.session} and \code{?is.session}
#'  for more details.
#' @param meanTol,varTol Tolerance for mean and variance of each data location. Locations which
#'  do not meet these thresholds are masked out of the analysis. Defaults: \code{1e-6}.
#' @param verbose Print messages counting how many locations are removed?
#'
#' @importFrom matrixStats colVars
#' @return A logical vector indicating valid vertices
#'
#' @export
make_mask <- function(data, meanTol=1e-6, varTol=1e-6, verbose=TRUE){

  # For each BOLD data matrix,
  mask_na <- mask_mean <- mask_var <- rep(TRUE, ncol(data[[1]]$BOLD))
  for (ss in seq(length(data))) {
    dss <- data[[ss]]$BOLD
    # Mark columns with any NA or NaN values for removal.
    dss_na <- is.na(dss) | is.nan(dss)
    mask_na[apply(dss_na, 2, any)] <- FALSE
    # Calculate means and variances of columns, except those with any NA or NaN.
    # Mark columns with mean/var falling under the thresholds for removal.
    mask_mean[mask_na][colMeans(dss[,mask_na,drop=FALSE]) < meanTol] <- FALSE
    mask_var[mask_na][matrixStats::colVars(dss[,mask_na,drop=FALSE]) < varTol] <- FALSE
  }

  # Print counts of locations removed, for each reason.
  if (verbose) {
    warn_part1 <- if (any(!mask_na)) { "additional locations" } else { "locations" }
    warn_part2 <- if (length(data) > 1) { " in at least one scan.\n" } else { ".\n" }
    if (any(!mask_na)) {
      cat("\t", sum(!mask_na), paste0("locations removed due to NA or NaN values", warn_part2))
    }
    # Do not include NA locations in count.
    mask_mean2 <- mask_mean | (!mask_na)
    if (any(!mask_mean2)) {
      cat("\t", sum(!mask_mean2), warn_part1, paste0("removed due to low mean", warn_part2))
    }
    # Do not include NA or low-mean locations in count.
    mask_var2 <- mask_var | (!mask_mean) | (!mask_na)
    if (any(!mask_var2)) {
      cat("\t", sum(!mask_var2), warn_part1, paste0("removed due to low variance", warn_part2))
    }
  }

  # Return composite mask.
  mask_na & mask_mean & mask_var
}