#' Retroactively mask locations from fit_bglm result.
#'
#' @param x The \code{"fit_bglm"} result
#' @param mask The mask to be applied to \code{x}. It's relative to the full
#'  mesh/voxel array, not \code{x}.
#
#' @return The masked result
#'
#' @keywords internal
retro_mask_fit_bglm <- function(x, mask){
  stopifnot(inherits(x, "fit_bglm"))
  spatial_type <- x$spatial$spatial_type

  mask_x <- switch(spatial_type,
    vertex = x$spatial$mask,
    voxel = x$spatial$labels != 0
  )

  nS <- length(x$session_names)
  nK <- length(x$field_names)
  nV_total <- length(mask_x)
  nV_input <- sum(mask_x)
  nT <- length(x$y) / nV_input / nS
  stopifnot(nT == round(nT))

  stopifnot(is.logical(mask))
  stopifnot(nV_total == length(mask))
  stopifnot(sum(mask) > 0)

  mask_new <- mask[mask_x]

  for (ss in seq(length(x$field_estimates))) {
    x$field_estimates[[ss]] <- x$field_estimates[[ss]][mask_new,,drop=FALSE]
  }

  for (ss in seq(length(x$RSS))) {
    x$RSS[[ss]] <- x$RSS[[ss]][mask_new]
  }

  if ("result_classical" %in% names(x)) {
    for (ss in seq(length(x$result_classical))) {
      x$result_classical[[ss]]$estimates <- x$result_classical[[ss]]$estimates[mask_new,,drop=FALSE]
      x$result_classical[[ss]]$SE_estimates <- x$result_classical[[ss]]$SE_estimates[mask_new,,drop=FALSE]
      x$result_classical[[ss]]$resids <- x$result_classical[[ss]]$resids[mask_new,,drop=FALSE]
      x$result_classical[[ss]]$RSS <- x$result_classical[[ss]]$RSS[mask_new]
    }
  }

  # Note: spde updated later

  # [TO DO] Ignore mask_qc right? It's for the user's interest; we don't use this.

  # Do this before `spatial` because the subcortex needs the old buffer mask.
  x$y <- c(matrix(x$y, ncol=nV_input)[,mask_new,drop=FALSE])
  for (ss in seq(length(x$X))) {
    print(dim(x$X[[ss]]))
    x$X[[ss]] <- if (spatial_type == "vertex") {
      x$X[[ss]][rep(mask_new, each=nT),rep(mask_new, each=nK),drop=FALSE]
    } else if (spatial_type == "voxel") {
      # something is wrong before here, when `mask_new` is not all TRUE.
      q <- x$spatial$buffer_mask;
      q[q] <- mask_new;
      q[!x$spatial$buffer_mask] <- TRUE;
      x$X[[ss]][rep(mask_new, each=nT),rep(q, each=nK),drop=FALSE]
    } else { stop() }
    print(dim(x$X[[ss]]))
  }

  if (spatial_type == "vertex") {
    x$spatial$mask <- mask
    x$spde$mesh <- retro_mask_mesh(x$spde$mesh, mask_new)
  } else if (spatial_type == "voxel") {
    x$spatial$labels[x$spatial$labels!=0][!mask_new] <- 0
    x$spatial$buffer_mask[x$spatial$buffer_mask] <- mask_new
    x$spatial$data_loc <- x$spatial$data_loc[mask_new]
    # [TO DO] in the future, leave labels alone?
  } else { stop() }

  cat("\n")

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
#' @param x The activations object
#' @param Masks The masks to be applied to each brain structure of \code{x}.
#' @return The masked result
#'
#' @keywords internal
retro_mask_act <- function(x, Masks){

  brainstructures <- names(Masks)
  nB <- length(Masks)
  nS <- length(x$activations[[1]])
  nG <- length(x$activations[[1]][[1]])

  for (bb in seq(nB)) {
    bs <- brainstructures[bb]
    # Get the mask to apply to the elements of `active`.
    spatial_type_bb <- x$spatial[[bb]]$spatial_type
    mask_bb <- switch(spatial_type_bb,
      surf = Masks[[bb]][x$spatial[[bb]]$mask],
      voxel = Masks[[bb]][x$spatial[[bb]]$labels[] != 0]
    )

    if (!all(mask_bb)) {
      message("\tRemoving ", sum(!mask_bb), " locations in ", bs, " model.")# for subject ", nn, ".")
    }

    # Apply the mask.
    for (ss in seq(nS)) {
      # Adjust `activations_xii`
      if (bs=="cortexL") {
        x$activations_xii[[ss]]$data$cortex_left <- x$activations_xii[[ss]]$data$cortex_left[mask_bb,,drop=FALSE]
        x$activations_xii[[ss]]$meta$cortex$medial_wall_mask$left[x$activations_xii[[ss]]$meta$cortex$medial_wall_mask$left] <- mask_bb
      } else if (bs=="cortexR") {
        x$activations_xii[[ss]]$data$cortex_right <- x$activations_xii[[ss]]$data$cortex_right[mask_bb,,drop=FALSE]
        x$activations_xii[[ss]]$meta$cortex$medial_wall_mask$right[x$activations_xii[[ss]]$meta$cortex$medial_wall_mask$right] <- mask_bb
      } else if (bs=="subcort") {
        x$activations_xii[[ss]]$data$subcort <- x$activations_xii[[ss]]$data$subcort[mask_bb,,drop=FALSE]
        x$activations_xii[[ss]]$meta$subcort$labels <- x$activations_xii[[ss]]$meta$subcort$labels[mask_bb]
        x$activations_xii[[ss]]$meta$subcort$mask[x$activations_xii[[ss]]$meta$subcort$mask] <- mask_bb
      } else { stop() }

      # Adjust `spatial`
      if (spatial_type_bb=="vertex") {
        x$spatial[[bb]]$mask <- mask_bb
      } else if (spatial_type_bb=="voxel") {
        # This will change...
        x$spatial[[bb]]$labels[x$spatial[[bb]]$labels != 0][!mask_bb] <- 0
      } else { stop() }

      # Adjust `activations`
      for (ss in seq(nS)) {
        for (gg in seq(nG)) {
          x$activations[[bb]][[ss]][[gg]]$active <- x$activations[[bb]][[ss]][[gg]]$active[mask_bb,,drop=FALSE]
        }
      }
    }
  }

  x
}
