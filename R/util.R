#' Plot BayesfMRI.spde objects
#'
#' @param object Object of class BayesfMRI.spde (see \code{help(create_spde_vol3D)})
#' @param colors (Optional) Vector of colors to represent each region.
#' @param alpha Transparency level.
#'
#' @return
#' @export
#' @importFrom geometry tetramesh
#' @import rgl
#' @import viridis viridis_pal
#'
#' @examples
plot.BayesfMRI.spde <- function(object, colors=NULL, alpha=0.5){
  if(class(object) != 'BayesfMRI.spde') stop('object argument must be a BayesfMRI.spde object. See help(create_spde_vol3D).')
  num_regions <- length(object$vertices)
  if(is.null(colors)) colors <- viridis_pal()(num_regions)
  if(length(colors) < num_regions) {
    colors <- rep(colors, length.out=num_regions)
    warning('Fewer colors specified than number of regions in the mesh.  Recycling colors to equal number of regions.')
  }
  for(ii in 1:num_regions){
    if(ii==1) tetramesh(object$faces[[ii]], object$vertices[[ii]], col=colors[ii], alpha=alpha)
    if(ii > 1) tetramesh(object$faces[[ii]], object$vertices[[ii]], clear=FALSE, col=colors[ii], alpha=alpha)

  }
}

#' Plot BayesGLM objects
#'
#' Summary method for class "BayesGLM"
#'
#' @param object an object of class "BayesGLM"
#' @param session_name NULL if BayesGLM object contains a single session; otherwise, the name of the session whose estimates to plot
#' @param pal If NULL, viridis palette with 64 colors will be used.  Otherwise, specify a vector of color names.
#' @export
#' @importFrom INLA plot.mesh
#' @import viridis
#' @method plot BayesGLM
plot.BayesGLM <- function(object, session_name=NULL, pal=NULL, ...)
{
  session_names <- names(result$beta_estimates)

  if((is.null(session_name)) & (length(session_names) > 1)) stop('If BayesGLM object includes multiple sessions, you must specify which session to plot.')
  if(!is.null(session_name) & !(session_name %in% session_names)) stop('I expect the session_names argument to be one of the session names of the BayesGLM object, but it is not.')

  if(is.null(session_name) & (length(session_names) == 1)) session_name <- session_names

  ind <- which(names(result$beta_estimates) == session_name) #which element of list
  est <- (result$beta_estimates)[[ind]]
  K <- ncol(est)


  if(is.null(pal)) {
    nColors <- 64
    pal <- viridis_pal(option='plasma', direction=1)(nColors)
  } else {
    if(min(areColors(pal)) < 1) stop('At least one of the elements of the pal argument is not a valid color representation.  See help(areColors).')
    nColors <- length(pal)
  }


  for(k in 1:K){
    x = est[,k]
    colindex <- as.integer(cut(x,breaks=nColors))

    #NEED TO CHECK WHICH TYPE OF BAYESGLM OBJECT (VOL OR CORTICAL) -- maybe use the mesh class?  or the spde_obj class?
    #plot(mesh_LH_s$mesh,col=pal[colindex], rgl=TRUE)

    tetramesh(object$spde_obj$faces, object$spde_obj$vertices, col=pal[colindex], clear=FALSE)

  }


}


#' Check whether each element of vector x is a valid color representation
#'
#' @param x Character vector
#'
#' @return A logical vector indicating which of the elements of x are valid color representations
#' @importFrom grDevices col2rgb
#' @export
#'
areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}


#' Summarise BayesGLM objects
#'
#' Summary method for class "BayesGLM"
#'
#' @param object an object of class "BayesGLM"
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary BayesGLM
summary.BayesGLM <- function(object, ...)
{
  out <- list()
  class(out) <- "summary.BayesGLM"
  out$sessions <- object$session_names
  out$betas <- object$beta_names
  out$call <- object$INLA_result$call
  out$inla.summary <- summary(object$model)
  return(out)
}


#' @param x an object of class "summary.BayesGLM"
#' @export
#' @method print summary.BayesGLM
#' @rdname summary.BayesGLM
print.summary.BayesGLM <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("Sessions: ", x$sessions,"\n")
  cat("Time used:\n", x$inla.summary$cpu.used)
}

#' @export
#' @method print BayesGLM
#' @rdname summary.BayesGLM
print.BayesGLM <- function(x, ...) {
  print.summary.BayesGLM(summary(x))
}

# TO DO: Add print and summary functions for session object (may need as.session function too, and move is.session here)



#' Find nonzero element in a matrix using 2-means clustering
#'
#' @param beta_est A vector or matrix of values from which values close to zero should be assigned a value of zero.
#'
#' @return A vector or matrix of the same dimension as beta_est in which values close to zero are assigned the value of zero. The closeness of a value to zero is found by performing two-means clustering on the absolute values of beta_est, and
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
#' @param x A vector consisting of all variables of interest for a single draw from a posterior distribution
#' @param b A scale parameter used to determine at what distance cluster centers are considered to be the same.
#'
#' @return The number of nonzero values detected within x
#' @export
#'
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

#' Array split
#'
#' Generalizes the `split` function to work with arrays of arbitrary dimension greater than 2.
#'
#' @param X An array of dimension greater than 1
#' @param margin The margin at which the array should be split
#'
#' @return A list of subarrays split along the `margin` dimension
#' @export
#'
asplit <- function(X,margin) {
  if(length(dim(X)) < 2) stop("This will only work matrices or arrays of dimension 2 or higher.")
  if(margin > length(dim(X))) stop("Dimension of X is less than margin.")
  X_dim <- dim(X)
  split_list <- split(X,slice.index(X,margin))
  output <- lapply(split_list,function(j){
    out <- array(j,dim = X_dim[-margin])
    return(out)
  })
  return(output)
}
