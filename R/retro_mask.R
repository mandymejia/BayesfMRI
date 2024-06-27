#' Retroactively mask locations from fit_bglm result.
#'
#' Work in progress.
#'
#' @param x The \code{"fit_bglm"} result
#' @param mask The mask to be applied to \code{x} (on top of any masks already
#'  applied to it.)
#' @return The masked result
#'
#' @keywords internal
retro_mask_BGLM <- function(x, mask){
  stopifnot(inherits(x, "fit_bglm"))
  nS <- length(x$session_names)
  nK <- length(x$field_names)
  nV <- sum(x$mask)
  nT <- length(x$y) / nV / nS
  stopifnot(nT == round(nT))

  stopifnot(is.logical(mask))
  stopifnot(nV == length(mask))
  stopifnot(sum(mask) > 0)

  mask2 <- x$mask
  x$mask[x$mask][!mask] <- FALSE
  x$mask[!mask2] <- FALSE

  for (ii in seq(length(x$field_estimates))) {
    x$field_estimates[[ii]][!mask2,] <- NA
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

#' Retroactively mask activations
#'
#' @param x The activation list
#' @param mask The mask to be applied to \code{x} (on top of any masks already
#'  applied to it.)
#' @return The masked result
#'
#' @keywords internal
retro_mask_act <- function(x){

  Masks <- intersect_mask_act(x)
  # [TO DO] tell the user about `Mask`, how many /bs
  brainstructures <- names(Masks)
  nB <- length(Masks)

  nN <- length(x)
  nS <- length(x[[1]]$activations[[1]])
  nG <- length(x[[1]]$activations[[1]][[1]])

  for (bb in seq(nB)) {
    bs <- brainstructures[bb]
    # Get the mask to apply to the elements of `active`.
    spatial_type_bb <- x[[1]]$spatial[[bb]]$spatial_type

    for (nn in seq(nN)) {
      mask_bb <- switch(spatial_type_bb,
        surf = Masks[[bb]][x[[nn]]$spatial[[bb]]$mask],
        voxel = Masks[[bb]][x[[nn]]$spatial[[bb]]$labels[] != 0]
      )

      if (!all(mask_bb)) {
        message("Intersection mask: removing ", sum(!mask_bb), " locations in ", bs, " model for subject ", nn, ".")
      }

      # Apply the mask.
      for (ss in seq(nS)) {
        # Adjust `activations_xii`
        if (bs=="cortexL") {
          x[[nn]]$activations_xii[[ss]]$data$cortex_left <- x[[nn]]$activations_xii[[ss]]$data$cortex_left[mask_bb,,drop=FALSE]
          x[[nn]]$activations_xii[[ss]]$meta$cortex$medial_wall_mask$left[x[[nn]]$activations_xii[[ss]]$meta$cortex$medial_wall_mask$left] <- mask_bb
        } else if (bs=="cortexR") {
          x[[nn]]$activations_xii[[ss]]$data$cortex_right <- x[[nn]]$activations_xii[[ss]]$data$cortex_right[mask_bb,,drop=FALSE]
          x[[nn]]$activations_xii[[ss]]$meta$cortex$medial_wall_mask$right[x[[nn]]$activations_xii[[ss]]$meta$cortex$medial_wall_mask$right] <- mask_bb
        } else if (bs=="subcort") {
          x[[nn]]$activations_xii[[ss]]$data$subcort <- x[[nn]]$activations_xii[[ss]]$data$subcort[mask_bb,,drop=FALSE]
          x[[nn]]$activations_xii[[ss]]$meta$subcort$labels <- x[[nn]]$activations_xii[[ss]]$meta$subcort$labels[mask_bb]
          x[[nn]]$activations_xii[[ss]]$meta$subcort$mask[x[[nn]]$activations_xii[[ss]]$meta$subcort$mask] <- mask_bb
        } else { stop() }

        # Adjust `spatial`
        if (spatial_type_bb=="surf") {
          x[[nn]]$spatial[[bb]]$mask <- mask_bb
        } else if (spatial_type_bb=="voxel") {
          x[[nn]]$spatial[[bb]]$labels <- x[[nn]]$spatial[[bb]]$labels[mask_bb]
          x[[nn]]$spatial[[bb]]$mask[x[[nn]]$spatial[[bb]]$mask] <- mask_bb
        } else { stop() }

        # Adjust `activations`
        for (ss in seq(nS)) {
          for (gg in seq(nG)) {
            x[[nn]]$activations[[bb]][[ss]][[gg]]$active <- x[[nn]]$activations[[bb]][[ss]][[gg]]$active[mask_bb,,drop=FALSE]
          }
        }
      }
    }
  }

  x
}
