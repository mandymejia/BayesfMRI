intersect_mask_act <- function(act_list) {
  brainstructures <- names(act_list[[1]]$spatial)
  nB <- length(brainstructures)

  # [TO DO] check that brainstructures match across sessions.

  # Get intersection mask for each brainstructure.
  Masks <- setNames(vector("list", nB), brainstructures)
  for (bb in seq(nB)) {
    bs <- brainstructures[bb]
    spatial_bb <- act_list[[1]]$spatial[[bs]]
    if("buffer_mask" %in% names(spatial_bb)) spatial_type_bb <- "voxel"
    if("surf" %in% names(spatial_bb)) spatial_type_bb <- "surf"
    if (spatial_type_bb=="surf") {
      masks <- do.call(rbind, lapply(act_list, function(x){
        x$spatial[[bs]]$mask
      }))

      # Identify non-empty locations with identical labels across sessions.
      Masks[[bb]] <- apply(masks, 2, all)

    } else if (spatial_type_bb=="voxel") {
      # Get vectorized, unmasked labels for every session.
      labels <- do.call(rbind, lapply(act_list, function(x){
        q <- x$spatial[[bs]]$mask * 1
        q[q==0] <- as.numeric(x$spatial[[bs]]$labels)
      }))

      # Identify non-empty locations with identical labels across sessions.
      Masks[[bb]] <- (colVars(labels)==0) & (labels[1,]!=0)
    } else { stop() }
  }

  Masks
}
