#' Retroactively mask locations from BayesGLM result.
#'
#' Work in progress.
#'
#' @param x The BayesGLM result
#' @param mask The mask to be applied to \code{x} (on top of any masks already
#'  applied to it.)
#' @return The masked result
#'
#' @keywords internal
retro_mask_BGLM <- function(x, mask){
  stopifnot(inherits(x, "BayesGLM"))
  nS <- length(x$session_names)
  nK <- length(x$task_names)
  nV <- sum(x$mask)
  nT <- length(x$y) / nV / nS
  stopifnot(nT == round(nT))

  stopifnot(is.logical(mask))
  stopifnot(nV == length(mask))
  stopifnot(sum(mask) > 0)

  mask2 <- x$mask
  x$mask[x$mask][!mask] <- FALSE
  x$mask[!mask2] <- FALSE

  for (ii in seq(length(x$task_estimates))) {
    x$task_estimates[[ii]][!mask2,] <- NA
  }

  if ("result_classical" %in% names(x)) {
    for (ii in seq(length(x$result_classical))) {
      x$result_classical[[ii]]$estimates[!mask2,] <- NA
      x$result_classical[[ii]]$SE_estimates[!mask2,] <- NA
      x$result_classical[[ii]]$resids <- x$result_classical[[ii]]$resids[mask,]
      x$result_classical[[ii]]$mask[!mask2] <- FALSE
    }
  }

  x$mesh <- retro_mask_mesh(x$mesh, mask)

  x$y <- c(matrix(x$y, ncol=nV)[,mask])
  for (ii in seq(length(x$X))) {
    x$X[[ii]] <- x$X[[ii]][rep(mask, each=nT),rep(mask, each=nK)]
  }

  x
}

#' Retroactively mask locations from mesh.
#'
#' Work in progress.
#'
#' @param x The mesh
#' @param mask The mask to be applied to \code{x} (on top of any masks already
#'  applied to it.)
#' @return The masked result
#'
#' @keywords internal
retro_mask_mesh <- function(x, mask){
  stopifnot(inherits(x, "inla.mesh"))
  stopifnot(is.logical(mask))

  # Determine which faces to keep.
  face_ok <- rowSums(matrix(x$graph$tv %in% which(mask), ncol=3)) == 3
  # Re-number the faces.
  fidx <- vector("numeric", length(mask))
  fidx[mask] <- seq(sum(mask))
  faces <- matrix(fidx[c(x$graph$tv[face_ok,])], ncol=3)
  # Make the new mesh.
  make_mesh(vertices=x$loc[mask,], faces=faces)
}
