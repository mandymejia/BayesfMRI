#' SPDE from mesh model
#'
#' @param spatial See \code{BayesGLM0}.
#' @return List
#' @keywords internal
SPDE_from_mesh <- function(spatial){
  surf <- spatial$surf
  mask <- spatial$mask
  nV <- nrow(surf$vertices)

  # Create INLA mesh.
  mesh_full <- make_mesh(surf$vertices, surf$faces)

  # [TO DO] change strategy; use a full mesh.
  # Create submesh of in-mask locations.
  # Update masks (sometimes `submesh` will exclude additional vertices).
  if (!all(mask)) {
    mesh <- submesh(mask, mesh_full)
    mask_new <- !is.na(mesh$idx$loc)
    mask_new_diff <- mask_new[mask]
    mask[mask] <- mask_new_diff
    mesh$idx$loc <- mesh$idx$loc[mask]
  }

  list(
    mesh_full = mesh_full,
    mesh = mesh,
    mask_new_diff = mask_new_diff,
    spde = INLA::inla.spde2.matern(mesh),
    spatial = list(surf=surf, mask=mask)
  )
}
