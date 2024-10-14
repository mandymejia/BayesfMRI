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

  maskIn_new <- mask[x$spatial$maskIn]
  maskMdat_new <- mask[x$spatial$maskMdat]

  if (any(!maskIn_new)) {
    for (ss in seq(nS)) {
      # Length nV_input of this session (not intersection)
      x$field_estimates[[ss]] <- x$field_estimates[[ss]][maskIn_new,,drop=FALSE]
      x$RSS[[ss]] <- x$RSS[[ss]][maskIn_new]
      if ("result_classical" %in% names(x)) {
        x$result_classical[[ss]]$estimates <- x$result_classical[[ss]]$estimates[maskIn_new,,drop=FALSE]
        x$result_classical[[ss]]$SE_estimates <- x$result_classical[[ss]]$SE_estimates[maskIn_new,,drop=FALSE]
        x$result_classical[[ss]]$resids <- x$result_classical[[ss]]$resids[maskIn_new,,drop=FALSE]
        x$result_classical[[ss]]$RSS <- x$result_classical[[ss]]$RSS[maskIn_new]
      }
      x$BOLD_QC$mask <- x$BOLD_QC$mask[maskIn_new]
      x$BOLD_QC$mask_na <- x$BOLD_QC$mask_na[maskIn_new]
      x$BOLD_QC$mask_mean <- x$BOLD_QC$mask_mean[maskIn_new]
      x$BOLD_QC$mask_var <- x$BOLD_QC$mask_var[maskIn_new]
      x$BOLD_QC$mask_snr <- x$BOLD_QC$mask_snr[maskIn_new]
    }
  }

  if (any(!maskMdat_new)) {
    for (ss in seq(nS)) {
      # Length nV_mdata of this session (not intersection)
      if ("prewhiten_info" %in% names(x)) {
        x$prewhiten_info$AR_coefs_avg <- x$prewhiten_info$AR_coefs_avg[maskMdat_new,,drop=FALSE]
        x$prewhiten_info$var_avg <- x$prewhiten_info$var_avg[maskMdat_new,drop=FALSE]
      }
    }

    # Notes for `BayesGLM2`:
    #   `spde` has been updated prior to `retro_mask_fit_bglm` call, and will
    #     be set after this function call.
    #   `spatial` will be determined later after this function call.

    # Do this before `spatial` because the subcortex needs the old buffer mask.
    x$y <- c(matrix(x$y, ncol=nV_input)[,maskMdat_new,drop=FALSE])
    for (ss in seq(length(x$X))) {
      x$X[[ss]] <- x$X[[ss]][rep(maskMdat_new, each=nT),,drop=FALSE]
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

  brainstructures <- names(Masks$Mdat)
  nB <- length(Masks$Mdat)

  one_bs_act <- if (inherits(x, "act_fit_bglm")) {
    x$activations
  } else if (inherits(x, "act_BGLM")) {
    x$activations[[1]]
  }
  nS <- length(one_bs_act)
  nG <- length(one_bs_act[[1]])

  for (bb in seq(nB)) {
    bs <- brainstructures[bb]
    # Get the mask to apply to the elements of `active`.
    spatial_type_bb <- x$spatial[[bb]]$spatial_type
    mask_bb <- Masks$Mdat[[bb]][as.logical(x$spatial[[bb]]$maskIn)]

    if (!all(mask_bb)) {
      message("\tRemoving ", sum(!mask_bb), " locations in ", bs, " model.")# for subject ", nn, ".")
    }

    # Apply the mask.
    for (ss in seq(nS)) {
      # Adjust `activations_xii`
      if (inherits(x, "act_BGLM")) {
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
      }

      # Adjust `activations`
      for (ss in seq(nS)) {
        for (gg in seq(nG)) {
          if (inherits(x, "act_BGLM")) {
            x$activations[[bb]][[ss]][[gg]]$active <- x$activations[[bb]][[ss]][[gg]]$active[mask_bb,,drop=FALSE]
          } else {
            x$activations[[ss]][[gg]]$active <- x$activations[[ss]][[gg]]$active[mask_bb,,drop=FALSE]
          }
        }
      }

      # Adjust `spatial`
      mask_bb <- Masks$Mdat[[bb]][as.logical(x$spatial[[bb]]$maskMdat)]
      x$spatial[[bb]]$maskMdat[x$spatial[[bb]]$maskMdat] <- mask_bb
      if (spatial_type_bb=="voxel") {
        x$spatial[[bb]]$labsMdat <- x$spatial[[bb]]$labsMdat[mask_bb]
      }
    }
  }

  x
}
