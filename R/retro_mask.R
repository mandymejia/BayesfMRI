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

  nS <- length(x$session_names)
  nK <- length(x$field_names)
  nV_total <- length(x$spatial$maskMdat)
  nV_input <- sum(x$spatial$maskMdat)
  nT <- length(x$y) / nV_input / nS
  stopifnot(nT == round(nT))

  stopifnot(is.logical(mask))
  stopifnot(nV_total == length(mask))
  stopifnot(sum(mask) > 0)

  mask_new <- mask[x$spatial$maskMdat]

  if (all(mask_new)) { return(x) }

  for (ss in seq(nS)) {
    x$field_estimates[[ss]] <- x$field_estimates[[ss]][mask_new,,drop=FALSE]
    x$RSS[[ss]] <- x$RSS[[ss]][mask_new]
    if ("result_classical" %in% names(x)) {
      x$result_classical[[ss]]$estimates <- x$result_classical[[ss]]$estimates[mask_new,,drop=FALSE]
      x$result_classical[[ss]]$SE_estimates <- x$result_classical[[ss]]$SE_estimates[mask_new,,drop=FALSE]
      x$result_classical[[ss]]$resids <- x$result_classical[[ss]]$resids[mask_new,,drop=FALSE]
      x$result_classical[[ss]]$RSS <- x$result_classical[[ss]]$RSS[mask_new]
    }
    x$BOLD_QC$mask <- x$BOLD_QC$mask[mask_new]
    x$BOLD_QC$mask_na <- x$BOLD_QC$mask_na[mask_new]
    x$BOLD_QC$mask_mean <- x$BOLD_QC$mask_mean[mask_new]
    x$BOLD_QC$mask_var <- x$BOLD_QC$mask_var[mask_new]
    x$BOLD_QC$mask_snr <- x$BOLD_QC$mask_snr[mask_new]
    if ("prewhiten_info" %in% names(x)) {
      x$prewhiten_info$AR_coefs_avg <- x$prewhiten_info$AR_coefs_avg[mask_new,,drop=FALSE]
      x$prewhiten_info$var_avg <- x$prewhiten_info$var_avg[mask_new,drop=FALSE]
    }
  }
  mask_new
  # Note: `spde` updated later

  # Do this before `spatial` because the subcortex needs the old buffer mask.
  x$y <- c(matrix(x$y, ncol=nV_input)[,mask_new,drop=FALSE])
  for (ss in seq(length(x$X))) {
    x$X[[ss]] <- x$X[[ss]][rep(mask_new, each=nT),,drop=FALSE]
    if (spatial_type == "voxel") {
      # Remove dropped locations from `X`.
      to_drop_Mdat <- which(!mask[x$spatial$maskMdat])
      if (length(to_drop_Mdat)!=0) {
        to_drop_spde <-x$spatial$Mmap[to_drop_Mdat]
        id_drop_spde <- seq(x$spde$n.spde) %in% to_drop_spde
        x$X[[ss]] <- x$X[[ss]][,!rep(id_drop_spde, times=nK),drop=FALSE]
      }
    }
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
    mask_bb <- Masks[[bb]][as.logical(x$spatial[[bb]]$maskMdat)]

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
      x$spatial[[bb]]$maskMdat[x$spatial[[bb]]$maskMdat] <- mask_bb
      if (spatial_type_bb=="voxel") {
        x$spatial[[bb]]$labsMdat <- x$spatial[[bb]]$labsMdat[mask_bb]
      }

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
