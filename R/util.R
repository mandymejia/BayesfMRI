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
#   out$betas <- object$beta_names
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



#' Find nonzero element in a matrix using 2-means clustering
#'
#' @param beta_est A vector or matrix of values from which values close to zero should be assigned a value of zero.
#'
#' @return A vector or matrix of the same dimension as beta_est in which values close to zero are assigned the value of zero. The closeness of a value to zero is found by performing two-means clustering on the absolute values of beta_est, and
#' @export
#' @importFrom stats kmeans
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
#' @importFrom stats quantile median
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

#' Boundary Mask
#'
#' Identify the vertices within `boundary_width` edges of the input mask. The
#'  faces must be triangular.
#'
#' @param faces a V x 3 matrix of integers. Each row defines a face by the index
#'  of three vertices.
#' @inheritParams mask_Param_vertices
#' @param boundary_width a positive integer. Vertices no more than this number
#'  of edges from any vertex in the input mask will be placed in the boundary mask.
#' @return The boundary mask, a length-V logical vector. TRUE indicates vertices
#'  within the boundary mask.
#'
boundary_mask <- function(faces, mask, boundary_width){
  s <- ncol(faces)
  v <- max(faces)
  # For quads, boundary_mask() would count opposite vertices on a face as
  #   adjacent--that's probably not desired.
  stopifnot(s == 3)

  stopifnot(boundary_width > 0)

  boundary_mask <- rep(FALSE, v)
  # Begin with the input mask.
  verts_adj_previous <- which(mask)
  for (ii in seq(1, boundary_width)) {
    # Identify vertices not in the mask, but adjacent to it.
    # Adjacency is defined by sharing a face.
    faces_nmask <- rowSums(matrix(faces %in% verts_adj_previous, ncol=s))
    faces_adj <- faces_nmask > 0 & faces_nmask < s
    verts_adj <- unique(as.vector(faces[faces_adj,]))
    verts_adj <- verts_adj[!(verts_adj %in% verts_adj_previous)]
    # Add those vertices to the boundary mask, and use them as the mask in
    #   the next iteration.
    boundary_mask[verts_adj] <- TRUE
    verts_adj_previous <- verts_adj
  }

  boundary_mask
}

#' Match user inputs to expected values
#'
#' Match each user input to an expected/allowed value. Raise a warning if either
#'  several user inputs match the same expected value, or at least one could not
#'  be matched to any expected value. \code{ciftiTools} uses this function to
#'  match keyword arguments for a function call. Another use is to match
#'  brainstructure labels ("left", "right", or "subcortical").
#'
#' @param user Character vector of user input. These will be matched to
#'  \code{expected} using \code{match.arg()}.
#' @param expected Character vector of expected/allowed values.
#' @param fail_action If any value in \code{user} could not be
#'  matched, or repeated matches occured, what should happen? Possible values
#'  are \code{"stop"} (default; raises an error), \code{"warning"}, and
#'  \code{"nothing"}.
#' @param user_value_label How to refer to the user input in a stop or warning
#'  message. If \code{NULL}, no label is used.
#'
#' @return The matched user inputs.
#'
#' @keywords internal
#'
match_input <- function(
  user, expected,
  fail_action=c("stop", "warning", "message", "nothing"),
  user_value_label=NULL) {

  fail_action <- match.arg(
    fail_action,
    c("stop", "warning", "message", "nothing")
  )
  unrecognized_FUN <- switch(fail_action,
                             stop=stop,
                             warning=warning,
                             message=message,
                             nothing=invisible
  )

  if (!is.null(user_value_label)) {
    user_value_label <- paste0("\"", user_value_label, "\" ")
  }
  msg <- paste0(
    "The user-input values ", user_value_label,
    "did not match their expected values. ",
    "Either several matched the same value, ",
    "or at least one did not match any.\n\n",
    "The user inputs were:\n",
    "\t\"", paste0(user, collapse="\", \""), "\".\n",
    "The expected values were:\n",
    "\t\"", paste0(expected, collapse="\", \""), "\".\n"
  )

  tryCatch(
    {
      matched <- match.arg(user, expected, several.ok=TRUE)
      if (length(matched) != length(user)) { stop() }
      return(matched)
    },
    error = function(e) {
      unrecognized_FUN(msg)
    },
    finally = {
    }
  )

  invisible(NULL)
}
